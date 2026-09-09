"""Asymptotic McDonald-Kreitman test implementation.

Based on Messer & Petrov (2013) PNAS.
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Callable

import numpy as np
from scipy import optimize

from mkado.analysis.statistics import omega_decomposition
from mkado.core.alignment import AlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


@dataclass
class PolymorphismData:
    """Intermediate data from a single gene for aggregation."""

    polymorphisms: list[tuple[float, str]] = field(default_factory=list)
    """List of (derived_frequency, type) where type is 'N' or 'S'."""
    dn: int = 0
    ds: int = 0
    gene_id: str = ""
    ln: float | None = None
    """Total non-synonymous sites (Nei-Gojobori) over analyzed codons."""
    ls: float | None = None
    """Total synonymous sites (Nei-Gojobori) over analyzed codons."""


def sum_site_totals(
    gene_data: list[PolymorphismData],
) -> tuple[float | None, float | None]:
    """Sum per-gene Nei-Gojobori (Ln, Ls); return ``(None, None)`` if any gene lacks them."""
    if not gene_data or any(g.ln is None or g.ls is None for g in gene_data):
        return (None, None)
    return (
        sum(g.ln for g in gene_data),  # type: ignore[misc]
        sum(g.ls for g in gene_data),  # type: ignore[misc]
    )


@dataclass
class AggregatedSFS:
    """Aggregated site frequency spectrum across genes."""

    bin_edges: np.ndarray = field(default_factory=lambda: np.array([]))
    pn_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    """Pn per bin (semantics depend on ``sfs_mode``: per-bin or cumulative tail)."""
    ps_counts: np.ndarray = field(default_factory=lambda: np.array([]))
    """Ps per bin (semantics depend on ``sfs_mode``: per-bin or cumulative tail)."""
    dn_total: int = 0
    ds_total: int = 0
    num_genes: int = 0
    ln_total: float | None = None
    """Sum of per-gene Nei-Gojobori non-synonymous sites."""
    ls_total: float | None = None
    """Sum of per-gene Nei-Gojobori synonymous sites."""
    pn_total: int = 0
    """Total non-synonymous polymorphisms (raw count, mode-invariant)."""
    ps_total: int = 0
    """Total synonymous polymorphisms (raw count, mode-invariant)."""
    sfs_mode: str = "at"
    """SFS construction mode used to build ``pn_counts``/``ps_counts``."""


@dataclass
class AsymptoticMKResult:
    """Results from an asymptotic McDonald-Kreitman test.

    ``Pn_total`` and ``Ps_total`` follow the standard MK definition: ingroup
    polymorphic sites whose derived allele causes a non-synonymous (Pn) or
    synonymous (Ps) substitution. The asymptotic test uses the full SFS
    with no ``min_frequency`` filter; frequency cutoffs apply only to the
    ``alpha(x)`` curve fit, not to the SFS counts themselves.
    """

    # Frequency bins and alpha estimates
    frequency_bins: list[float] = field(default_factory=list)
    """All frequency bin centers."""
    alpha_by_freq: list[float] = field(default_factory=list)
    """Alpha values for bins with valid data (same length as alpha_x_values)."""
    alpha_x_values: list[float] = field(default_factory=list)
    """X values (frequency bin centers) corresponding to alpha_by_freq."""

    # Asymptotic alpha estimate (extrapolated to x=1)
    alpha_asymptotic: float = 0.0

    # 95% confidence interval
    ci_low: float = 0.0
    ci_high: float = 0.0

    # Fit parameters (alpha = a + b*exp(-c*x) for exponential, a + b*x for linear)
    fit_a: float = 0.0
    fit_b: float = 0.0
    fit_c: float = 0.0

    # Counts used in calculation
    dn: int = 0
    ds: int = 0

    # Aggregation info (0 = single gene, >0 = aggregated)
    num_genes: int = 0

    # Model type ("exponential" or "linear")
    model_type: str = "exponential"

    pn_total: int = 0
    """Total non-synonymous polymorphisms (sum across all SFS bins)."""
    ps_total: int = 0
    """Total synonymous polymorphisms (sum across all SFS bins)."""
    ln: float | None = None
    """Total non-synonymous sites (Nei-Gojobori)."""
    ls: float | None = None
    """Total synonymous sites (Nei-Gojobori)."""
    omega: float | None = None
    """``(Dn/Ds) * (Ls/Ln)`` — dN/dS ratio."""
    omega_a: float | None = None
    """Adaptive component ``alpha_asymptotic * omega``."""
    omega_na: float | None = None
    """Non-adaptive component ``(1 - alpha_asymptotic) * omega``."""
    ci_method: str = "monte-carlo"
    """Method used to compute ``ci_low``/``ci_high``: ``"monte-carlo"`` samples
    parameters from the curve-fit covariance matrix; ``"bootstrap"`` resamples
    polymorphisms with replacement and refits per replicate."""
    sfs_mode: str = "at"
    """SFS construction mode: ``"at"`` (Messer & Petrov 2013, count per bin) or
    ``"above"`` (Uricchio et al. 2019, inclusive right-tail cumulative count)."""

    # Dn, Ds, Ln, Ls are constants under the asymptotic Monte Carlo procedure,
    # so omega has no sampling distribution. The CIs on omega_a and omega_na
    # are pure arithmetic transforms of the existing alpha CI: no need to
    # store them.
    @property
    def omega_a_ci_low(self) -> float | None:
        """Lower 95% CI for ``omega_a`` = ``ci_low * omega``."""
        if self.omega is None:
            return None
        return self.ci_low * self.omega

    @property
    def omega_a_ci_high(self) -> float | None:
        """Upper 95% CI for ``omega_a`` = ``ci_high * omega``."""
        if self.omega is None:
            return None
        return self.ci_high * self.omega

    @property
    def omega_na_ci_low(self) -> float | None:
        """Lower 95% CI for ``omega_na`` = ``(1 - ci_high) * omega``.

        Percentiles flip because ``(1 - alpha)`` is monotonically decreasing.
        """
        if self.omega is None:
            return None
        return (1.0 - self.ci_high) * self.omega

    @property
    def omega_na_ci_high(self) -> float | None:
        """Upper 95% CI for ``omega_na`` = ``(1 - ci_low) * omega``."""
        if self.omega is None:
            return None
        return (1.0 - self.ci_low) * self.omega

    def __str__(self) -> str:
        lines = [
            "Asymptotic MK Test Results:",
            f"  Asymptotic α: {self.alpha_asymptotic:.4f} "
            f"(95% CI [{self.ci_method}]: {self.ci_low:.4f} - {self.ci_high:.4f})",
            f"  Divergence: Dn={self.dn}, Ds={self.ds}",
            f"  SFS mode: {self.sfs_mode}",
        ]
        if self.pn_total > 0 or self.ps_total > 0:
            lines.append(f"  Polymorphism: Pn={self.pn_total}, Ps={self.ps_total}")
        if self.ln is not None and self.ls is not None:
            lines.append(f"  Sites: Ln={self.ln:.2f}, Ls={self.ls:.2f}")
        if self.omega is not None:
            omega_a_str = f"{self.omega_a:.4f}" if self.omega_a is not None else "NA"
            omega_na_str = f"{self.omega_na:.4f}" if self.omega_na is not None else "NA"
            lines.append(
                f"  omega: {self.omega:.4f} (omega_a={omega_a_str}, omega_na={omega_na_str})"
            )
            if self.omega_a_ci_low is not None and self.omega_a_ci_high is not None:
                lines.append(
                    f"    omega_a 95% CI:  ({self.omega_a_ci_low:.4f}, {self.omega_a_ci_high:.4f})"
                )
            if self.omega_na_ci_low is not None and self.omega_na_ci_high is not None:
                lines.append(
                    f"    omega_na 95% CI: ({self.omega_na_ci_low:.4f}, "
                    f"{self.omega_na_ci_high:.4f})"
                )
        if self.num_genes > 0:
            lines.append(f"  Genes aggregated: {self.num_genes}")
        if self.model_type == "exponential":
            lines.append(
                f"  Fit ({self.model_type}): α(x) = {self.fit_a:.4f} + "
                f"({self.fit_b:.4f}) * exp(-{self.fit_c:.4f} * x)"
            )
        else:
            lines.append(
                f"  Fit ({self.model_type}): α(x) = {self.fit_a:.4f} + {self.fit_b:.4f} * x"
            )
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        result = {
            "alpha_asymptotic": self.alpha_asymptotic,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "dn": self.dn,
            "ds": self.ds,
            "model_type": self.model_type,
            "fit_parameters": {
                "a": self.fit_a,
                "b": self.fit_b,
            },
            "frequency_bins": self.frequency_bins,
            "alpha_by_freq": self.alpha_by_freq,
            "alpha_x_values": self.alpha_x_values,
        }
        if self.model_type == "exponential":
            result["fit_parameters"]["c"] = self.fit_c
        if self.num_genes > 0:
            result["num_genes"] = self.num_genes
        if self.pn_total > 0 or self.ps_total > 0:
            result["pn_total"] = self.pn_total
            result["ps_total"] = self.ps_total
        result["ln"] = self.ln
        result["ls"] = self.ls
        result["omega"] = self.omega
        result["omega_a"] = self.omega_a
        result["omega_na"] = self.omega_na
        result["omega_a_ci_low"] = self.omega_a_ci_low
        result["omega_a_ci_high"] = self.omega_a_ci_high
        result["omega_na_ci_low"] = self.omega_na_ci_low
        result["omega_na_ci_high"] = self.omega_na_ci_high
        result["ci_method"] = self.ci_method
        result["sfs_mode"] = self.sfs_mode
        return result


def _attach_omega(
    result: AsymptoticMKResult,
    ln: float | None,
    ls: float | None,
) -> AsymptoticMKResult:
    """Populate ln/ls/omega point estimates on an asymptotic result in-place."""
    omega, omega_a, omega_na = omega_decomposition(
        result.dn, result.ds, ln, ls, result.alpha_asymptotic
    )
    result.ln = ln
    result.ls = ls
    result.omega = omega
    result.omega_a = omega_a
    result.omega_na = omega_na
    return result


_SFS_MODES: tuple[str, ...] = ("at", "above")


def _validate_sfs_mode(sfs_mode: str) -> None:
    """Reject any value other than 'at' or 'above'."""
    if sfs_mode not in _SFS_MODES:
        raise ValueError(f"sfs_mode must be one of {_SFS_MODES!r}, got {sfs_mode!r}")


def _apply_sfs_mode(pn: np.ndarray, ps: np.ndarray, sfs_mode: str) -> tuple[np.ndarray, np.ndarray]:
    """Transform per-bin counts into the requested SFS form.

    Under ``"at"`` mode the inputs are returned unchanged. Under ``"above"``
    mode each is replaced by its inclusive right-tail cumulative sum so that
    bin ``i`` holds the count in bins ``[i, end]`` (Uricchio et al. 2019).
    """
    if sfs_mode == "above":
        return np.cumsum(pn[::-1])[::-1], np.cumsum(ps[::-1])[::-1]
    return pn, ps


def _frequency_bin_edges(num_bins: int) -> np.ndarray:
    """Bin edges 0, 1/num_bins, ..., 1, each the nearest double to k/num_bins.

    ``np.linspace(0, 1, num_bins + 1)`` does not have this property for
    several ``num_bins``/``k`` combinations (e.g. ``num_bins=20``, edge 3
    comes back as 0.15000000000000002), which pushes a polymorphism at
    exactly that frequency into the wrong bin.
    """
    return np.arange(num_bins + 1) / num_bins


def _exponential_model(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Exponential model: α(x) = a + b * exp(-c * x)

    Following the asymptoticMK R package convention.
    For positive selection with deleterious mutations, b is typically negative
    (alpha increases with frequency).
    """
    return a + b * np.exp(-c * x)


def _exponential_model_jac(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    """Closed-form Jacobian of ``_exponential_model`` w.r.t. (a, b, c).

    Replaces ``scipy.optimize._numdiff.approx_derivative`` (~21% of bootstrap-CI
    wall time on the human dataset) with the analytic columns of shape (m, 3):
    ``[1, exp(-cx), -b·x·exp(-cx)]``.
    """
    x = np.asarray(x, dtype=float)
    e = np.exp(-c * x)
    out = np.empty((x.size, 3), dtype=float)
    out[:, 0] = 1.0
    out[:, 1] = e
    out[:, 2] = -b * x * e
    return out


def _linear_model(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Linear model: α(x) = a + b * x"""
    return a + b * x


def _linear_model_jac(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Closed-form Jacobian of ``_linear_model`` w.r.t. (a, b): columns [1, x]."""
    x = np.asarray(x, dtype=float)
    out = np.empty((x.size, 2), dtype=float)
    out[:, 0] = 1.0
    out[:, 1] = x
    return out


def _compute_ci_monte_carlo(
    popt: np.ndarray,
    pcov: np.ndarray,
    model_func: Callable[..., np.ndarray],
    n_sim: int = 10000,
    seed: int = 42,
) -> tuple[float, float]:
    """Compute CI via Monte Carlo simulation from covariance matrix.

    Args:
        popt: Optimal parameter values from curve_fit
        pcov: Covariance matrix from curve_fit
        model_func: Model function (_exponential_model or _linear_model)
        n_sim: Number of Monte Carlo simulations
        seed: Random seed for reproducibility

    Returns:
        Tuple of (ci_low, ci_high) for alpha at x=1
    """
    rng = np.random.default_rng(seed)

    # Handle case where covariance contains inf/nan
    if not np.all(np.isfinite(pcov)):
        alpha_at_1 = float(model_func(np.array([1.0]), *popt)[0])
        return (alpha_at_1, alpha_at_1)

    # Sample parameters from multivariate normal
    try:
        param_samples = rng.multivariate_normal(popt, pcov, size=n_sim)
    except np.linalg.LinAlgError:
        alpha_at_1 = float(model_func(np.array([1.0]), *popt)[0])
        return (alpha_at_1, alpha_at_1)

    # Evaluate model at x=1 for each sample
    alpha_samples = []
    for params in param_samples:
        try:
            alpha = float(model_func(np.array([1.0]), *params)[0])
            if np.isfinite(alpha):
                alpha_samples.append(alpha)
        except (ValueError, RuntimeWarning):
            continue

    if len(alpha_samples) < 100:
        alpha_at_1 = float(model_func(np.array([1.0]), *popt)[0])
        return (alpha_at_1, alpha_at_1)

    ci_low = float(np.percentile(alpha_samples, 2.5))
    ci_high = float(np.percentile(alpha_samples, 97.5))
    return (ci_low, ci_high)


def _bootstrap_one_replicate(
    sample: np.ndarray,
    bin_idx: np.ndarray,
    is_n: np.ndarray,
    num_bins: int,
    sfs_mode: str,
    ds: int,
    dn: int,
    bin_in_window: np.ndarray,
    bin_centers: np.ndarray,
    model_func: Callable[..., np.ndarray],
    p0: list[float],
    bounds: tuple[list[float], list[float]],
    model_jac: Callable[..., np.ndarray] | None = None,
) -> float | None:
    """Resample, rebin, refit α(x), evaluate at x=1. Returns None on failure."""
    sample_bins = bin_idx[sample]
    sample_n = is_n[sample]
    boot_pn = np.bincount(sample_bins[sample_n], minlength=num_bins).astype(float)
    boot_ps = np.bincount(sample_bins[~sample_n], minlength=num_bins).astype(float)
    boot_pn, boot_ps = _apply_sfs_mode(boot_pn, boot_ps, sfs_mode)

    # Per-bin alpha where boot_ps > 0 and inside the fitting window, falling
    # back to the full range if cutoffs are too restrictive (mirrors
    # _mask_to_frequency_cutoffs, the point-estimate path's equivalent).
    ps_positive = boot_ps > 0
    valid = _select_mask_with_fallback(ps_positive & bin_in_window, ps_positive)
    if valid.sum() < 3:
        return None
    ratio = (ds / dn) * (boot_pn[valid] / boot_ps[valid])
    boot_alphas = 1.0 - ratio
    boot_x = bin_centers[valid]

    try:
        popt_boot, _ = optimize.curve_fit(
            model_func,
            boot_x,
            boot_alphas,
            p0=p0,
            bounds=bounds,
            maxfev=5000,
            jac=model_jac,
        )
        alpha_at_1 = float(model_func(np.array([1.0]), *popt_boot)[0])
        if np.isfinite(alpha_at_1):
            return alpha_at_1
    except (RuntimeError, ValueError):
        pass
    return None


def _bootstrap_batch(args: tuple) -> list[float | None]:
    """Process a batch of bootstrap samples in one worker, returning a result list.

    Pickling shared state once per batch (rather than once per sample) is the
    point of batching: with ``len(samples_chunk) ≈ n_replicates / workers``,
    the per-task pickle cost is amortized across many replicates.
    """
    (
        samples_chunk,
        bin_idx,
        is_n,
        num_bins,
        sfs_mode,
        ds,
        dn,
        bin_in_window,
        bin_centers,
        model_func,
        p0,
        bounds,
        model_jac,
    ) = args
    return [
        _bootstrap_one_replicate(
            s,
            bin_idx,
            is_n,
            num_bins,
            sfs_mode,
            ds,
            dn,
            bin_in_window,
            bin_centers,
            model_func,
            p0,
            bounds,
            model_jac,
        )
        for s in samples_chunk
    ]


def _compute_ci_bootstrap(
    pooled_polymorphisms: list[tuple[float, str]],
    dn: int,
    ds: int,
    bin_edges: np.ndarray,
    bin_centers: np.ndarray,
    num_bins: int,
    model_func: Callable[..., np.ndarray],
    p0: list[float],
    bounds: tuple[list[float], list[float]],
    frequency_cutoffs: tuple[float, float],
    point_estimate: float,
    n_replicates: int = 100,
    seed: int = 42,
    sfs_mode: str = "at",
    workers: int = 1,
    model_jac: Callable[..., np.ndarray] | None = None,
) -> tuple[float, float]:
    """Compute CI by case-resampling the pooled polymorphism list.

    Each replicate draws ``len(pooled_polymorphisms)`` indices uniformly with
    replacement, rebins, recomputes per-bin α(x), refits ``model_func`` on the
    same model the point estimate selected, and evaluates at x=1. Replicates
    where the fit fails or fewer than 3 valid bins remain are dropped. If
    fewer than half of the replicates succeed, the CI degenerates to the
    point estimate.

    With ``workers > 1`` the per-replicate work is dispatched to a
    ``ProcessPoolExecutor`` in batches of ``n_replicates // workers``.
    Threads were tried first but ``scipy.optimize.curve_fit`` runs the
    Levenberg-Marquardt iteration in Python (releasing the GIL only inside
    individual numpy ops), so threading gives no speedup on these small
    fits. Processes pickle once per batch — fine because
    ``n_replicates >> workers``. Determinism is preserved: sample arrays
    are generated serially before dispatch, and ``np.percentile`` is
    order-invariant.
    """
    if workers < 1:
        raise ValueError(f"workers must be >= 1, got {workers}")

    n = len(pooled_polymorphisms)
    if n == 0 or dn == 0 or ds == 0 or n_replicates <= 0:
        return (point_estimate, point_estimate)

    # Vectorize: precompute frequency and is-N arrays once.
    freqs = np.fromiter((f for f, _ in pooled_polymorphisms), dtype=float, count=n)
    is_n = np.fromiter((t == "N" for _, t in pooled_polymorphisms), dtype=bool, count=n)

    # Precompute per-polymorphism bin indices once.
    bin_idx = np.searchsorted(bin_edges[1:], freqs, side="right")
    np.clip(bin_idx, 0, num_bins - 1, out=bin_idx)

    rng = np.random.default_rng(seed)
    low_cut, high_cut = frequency_cutoffs
    bin_in_window = (bin_centers >= low_cut) & (bin_centers <= high_cut)

    # Pre-generate sample arrays serially so the RNG is consumed deterministically;
    # parallel dispatch must not depend on worker scheduling for reproducibility.
    samples = [rng.integers(0, n, size=n) for _ in range(n_replicates)]

    if workers == 1:
        replicate_results = [
            _bootstrap_one_replicate(
                s,
                bin_idx,
                is_n,
                num_bins,
                sfs_mode,
                ds,
                dn,
                bin_in_window,
                bin_centers,
                model_func,
                p0,
                bounds,
                model_jac,
            )
            for s in samples
        ]
    else:
        # Round-robin chunks so each worker gets a similar mix of samples;
        # the SET of returned alphas is invariant to chunking, so np.percentile
        # gives byte-identical CI bounds regardless of worker count.
        chunks = [samples[i::workers] for i in range(workers)]
        batch_args = [
            (
                chunk,
                bin_idx,
                is_n,
                num_bins,
                sfs_mode,
                ds,
                dn,
                bin_in_window,
                bin_centers,
                model_func,
                p0,
                bounds,
                model_jac,
            )
            for chunk in chunks
        ]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            chunk_results = list(pool.map(_bootstrap_batch, batch_args))
        replicate_results = [r for chunk in chunk_results for r in chunk]

    bootstrap_alphas: list[float] = [r for r in replicate_results if r is not None]

    if len(bootstrap_alphas) < n_replicates // 2:
        logger.debug(
            "Bootstrap CI: only %d/%d replicates succeeded; falling back to point estimate",
            len(bootstrap_alphas),
            n_replicates,
        )
        return (point_estimate, point_estimate)

    return (
        float(np.percentile(bootstrap_alphas, 2.5)),
        float(np.percentile(bootstrap_alphas, 97.5)),
    )


def _compute_aic(n: int, rss: float, k: int) -> float:
    """Compute Akaike Information Criterion.

    Args:
        n: Number of data points
        rss: Residual sum of squares
        k: Number of parameters

    Returns:
        AIC value
    """
    if rss <= 0 or n <= k:
        return float("inf")
    return n * np.log(rss / n) + 2 * k


def _select_mask_with_fallback(
    candidate: np.ndarray, fallback: np.ndarray, min_points: int = 3
) -> np.ndarray:
    """Use `candidate` if it selects at least `min_points`, else `fallback`.

    Both the point-estimate and bootstrap fits trim to the requested
    frequency window when there's enough support to constrain a
    3-parameter exponential fit, and silently widen back to the full range
    otherwise, so a narrow --freq-cutoffs window doesn't starve the fit.
    """
    return candidate if candidate.sum() >= min_points else fallback


def _mask_to_frequency_cutoffs(
    x_data: np.ndarray, y_data: np.ndarray, frequency_cutoffs: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """Restrict (x_data, y_data) to frequency_cutoffs for curve fitting.

    Falls back to the full, untrimmed arrays when the cutoffs leave fewer
    than 3 points, since a 3-parameter exponential fit needs at least that
    many to be constrained. Shared by the aggregated and per-gene fitters
    so both trim bins the same way, whether the input is many genes pooled
    or a single gene.
    """
    low_cut, high_cut = frequency_cutoffs
    window_mask = (x_data >= low_cut) & (x_data <= high_cut)
    mask = _select_mask_with_fallback(window_mask, np.ones_like(window_mask))
    return x_data[mask], y_data[mask]


def extract_polymorphism_data(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    pool_polymorphisms: bool = False,
    gene_id: str = "",
    min_frequency: float = 0.0,
    no_singletons: bool = False,
) -> PolymorphismData:
    """Extract polymorphism and divergence data without curve fitting.

    This function extracts the raw data needed for asymptotic MK analysis,
    allowing multiple genes to be processed and aggregated before fitting.

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup: SequenceSet or path to FASTA file for outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        pool_polymorphisms: If True, consider sites polymorphic in either
            population (libsequence convention). Frequencies are still
            calculated from ingroup only.
        gene_id: Identifier for this gene
        min_frequency: Minimum derived allele frequency for polymorphisms
            (polymorphisms below this threshold are excluded)

    Returns:
        PolymorphismData containing polymorphisms with frequencies and divergence counts
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Create aligned pair
    pair = AlignedPair(ingroup=ingroup, outgroup=outgroup, genetic_code=code)

    # Count divergence
    dn = 0
    ds = 0

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            result = pair.classify_fixed_difference(codon_idx)
            if result is not None:
                nonsyn, syn = result
                dn += nonsyn
                ds += syn

    # Collect polymorphisms with their frequencies
    poly_data: list[tuple[float, str]] = []

    # Get polymorphic sites based on mode
    if pool_polymorphisms:
        poly_sites = pair.polymorphic_sites_pooled()
    else:
        poly_sites = pair.polymorphic_sites_ingroup()

    for codon_idx in poly_sites:
        codons = list(ingroup.codon_set_clean(codon_idx))
        if len(codons) < 2:
            continue

        state = ingroup.derived_state(outgroup, codon_idx)
        if state is None:
            continue
        derived_freq, derived_count = state

        if derived_freq <= 0 or derived_freq >= 1:
            continue

        if no_singletons and derived_count == 1:
            continue
        if derived_freq < min_frequency:
            continue

        # Classify the polymorphism
        result = pair.classify_polymorphism(codon_idx)
        if result is not None:
            nonsyn, syn = result
            for _ in range(nonsyn):
                poly_data.append((derived_freq, "N"))
            for _ in range(syn):
                poly_data.append((derived_freq, "S"))

    ln, ls = pair.count_total_sites()

    return PolymorphismData(
        polymorphisms=poly_data,
        dn=dn,
        ds=ds,
        gene_id=gene_id,
        ln=ln,
        ls=ls,
    )


def aggregate_polymorphism_data(
    gene_data: list[PolymorphismData],
    num_bins: int = 20,
    sfs_mode: str = "at",
) -> AggregatedSFS:
    """Aggregate polymorphism data from multiple genes into SFS bins.

    Args:
        gene_data: List of PolymorphismData from individual genes
        num_bins: Number of frequency bins
        sfs_mode: ``"at"`` (default, Messer & Petrov 2013) returns per-bin
            counts. ``"above"`` (Uricchio et al. 2019) returns the inclusive
            right-tail cumulative SFS — bin ``i`` holds the count of
            polymorphisms with frequency ``>= bin_edges[i]``.

    Returns:
        AggregatedSFS with per-bin Pn/Ps counts and total divergence
    """
    _validate_sfs_mode(sfs_mode)

    # Sum divergence across all genes
    dn_total = sum(g.dn for g in gene_data)
    ds_total = sum(g.ds for g in gene_data)
    ln_total, ls_total = sum_site_totals(gene_data)

    # Create frequency bins
    bin_edges = _frequency_bin_edges(num_bins)

    # Vectorize: collect all (freq, type) pairs into two parallel arrays, then
    # compute bin indices in one searchsorted and aggregate counts via bincount.
    # The previous per-polymorphism loop spent ~25% of aggregate_polymorphism_data
    # time inside a Python-level np.searchsorted call per item.
    n_total = sum(len(g.polymorphisms) for g in gene_data)
    if n_total == 0:
        pn_counts = np.zeros(num_bins)
        ps_counts = np.zeros(num_bins)
    else:
        freqs = np.empty(n_total, dtype=float)
        is_n = np.empty(n_total, dtype=bool)
        i = 0
        for g in gene_data:
            for freq, ptype in g.polymorphisms:
                freqs[i] = freq
                is_n[i] = ptype == "N"
                i += 1
        bin_idx = np.searchsorted(bin_edges[1:], freqs, side="right")
        np.clip(bin_idx, 0, num_bins - 1, out=bin_idx)
        pn_counts = np.bincount(bin_idx[is_n], minlength=num_bins).astype(float)
        ps_counts = np.bincount(bin_idx[~is_n], minlength=num_bins).astype(float)

    pn_total = int(pn_counts.sum())
    ps_total = int(ps_counts.sum())

    pn_counts, ps_counts = _apply_sfs_mode(pn_counts, ps_counts, sfs_mode)

    return AggregatedSFS(
        bin_edges=bin_edges,
        pn_counts=pn_counts,
        ps_counts=ps_counts,
        dn_total=dn_total,
        ds_total=ds_total,
        num_genes=len(gene_data),
        ln_total=ln_total,
        ls_total=ls_total,
        pn_total=pn_total,
        ps_total=ps_total,
        sfs_mode=sfs_mode,
    )


def asymptotic_mk_test_aggregated(
    gene_data: list[PolymorphismData],
    num_bins: int = 20,
    ci_replicates: int = 10000,
    frequency_cutoffs: tuple[float, float] = (0.1, 0.9),
    ci_method: str = "monte-carlo",
    seed: int = 42,
    sfs_mode: str = "at",
    workers: int = 1,
) -> AsymptoticMKResult:
    """Genome-wide asymptotic MK test on aggregated data.

    This follows the approach of Messer & Petrov (2013) and the asymptoticMK
    R package: aggregate polymorphism SFS data and divergence counts from
    many genes, then fit a single exponential curve to the aggregated data.

    Args:
        gene_data: List of PolymorphismData from individual genes
        num_bins: Number of frequency bins (default 20)
        ci_replicates: Number of CI replicates. For ``ci_method="monte-carlo"``
            (default), parameters are sampled from the curve-fit covariance —
            pass ~10000. For ``ci_method="bootstrap"``, each replicate refits
            the curve so 100-500 is typical.
        frequency_cutoffs: (low, high) frequency range for fitting (default 0.1-0.9)
        ci_method: ``"monte-carlo"`` (default, parametric MVN sampling from the
            curve-fit covariance) or ``"bootstrap"`` (case-resampling of the
            pooled polymorphism list with replacement, refit per replicate).
        seed: Random seed for reproducibility (default 42).

    Returns:
        AsymptoticMKResult with aggregated analysis
    """
    if ci_method not in ("monte-carlo", "bootstrap"):
        raise ValueError(f"ci_method must be 'monte-carlo' or 'bootstrap', got {ci_method!r}")
    _validate_sfs_mode(sfs_mode)

    # Aggregate the data
    agg = aggregate_polymorphism_data(gene_data, num_bins, sfs_mode=sfs_mode)

    dn = agg.dn_total
    ds = agg.ds_total
    bin_edges = agg.bin_edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Total polymorphisms (raw, mode-invariant — used for the simple-α fallback
    # and for output reporting, never for the per-bin α(x) curve).
    pn_total = agg.pn_total
    ps_total = agg.ps_total

    # Calculate per-bin alpha for each frequency bin (matching asymptoticMK R package
    # under "at" mode; under "above" mode this is α(>x) on the cumulative SFS).
    alpha_values = []
    valid_centers = []

    for i in range(num_bins):
        pn_bin = float(agg.pn_counts[i])
        ps_bin = float(agg.ps_counts[i])

        if ds > 0 and dn > 0 and ps_bin > 0:
            alpha_x = 1.0 - (ds / dn) * (pn_bin / ps_bin)
            alpha_values.append(alpha_x)
            valid_centers.append(bin_centers[i])

    # Not enough data for curve fitting
    if len(valid_centers) < 3:
        simple_alpha = None
        if ds > 0 and dn > 0 and ps_total > 0:
            simple_alpha = 1.0 - (ds * pn_total) / (dn * ps_total)

        result = AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_x_values=valid_centers,
            alpha_asymptotic=simple_alpha or 0.0,
            ci_low=simple_alpha or 0.0,
            ci_high=simple_alpha or 0.0,
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            pn_total=pn_total,
            ps_total=ps_total,
            ci_method=ci_method,
            sfs_mode=sfs_mode,
        )
        return _attach_omega(result, agg.ln_total, agg.ls_total)

    x_data = np.array(valid_centers)
    y_data = np.array(alpha_values)
    x_fit, y_fit = _mask_to_frequency_cutoffs(x_data, y_data, frequency_cutoffs)

    # Try exponential fit first
    exp_success = False
    exp_popt = None
    exp_pcov = None
    exp_alpha = 0.0
    exp_ci_width = float("inf")

    try:
        # Initial guesses: a is asymptotic value, b is negative for increasing alpha
        p0_exp = [y_fit[-1], y_fit[0] - y_fit[-1], 5.0]
        # b can be negative (typically is for increasing alpha with frequency)
        bounds_exp = ([-2.0, -2.0, 0.001], [2.0, 2.0, 100.0])

        exp_popt, exp_pcov = optimize.curve_fit(
            _exponential_model,
            x_fit,
            y_fit,
            p0=p0_exp,
            bounds=bounds_exp,
            maxfev=10000,
            jac=_exponential_model_jac,
        )

        a, b, c = exp_popt
        # α(x=1) = a + b * exp(-c)
        exp_alpha = float(a + b * np.exp(-c))

        # Compute CI via Monte Carlo
        ci_low_exp, ci_high_exp = _compute_ci_monte_carlo(
            exp_popt, exp_pcov, _exponential_model, ci_replicates
        )
        exp_ci_width = ci_high_exp - ci_low_exp
        exp_success = True

    except (RuntimeError, ValueError) as exc:
        logger.debug("Exponential model fit failed: %s", exc)

    # Try linear fit
    lin_success = False
    lin_popt = None
    lin_pcov = None
    lin_alpha = 0.0

    try:
        p0_lin = [y_fit[0], y_fit[-1] - y_fit[0]]
        bounds_lin = ([-2.0, -2.0], [2.0, 2.0])

        lin_popt, lin_pcov = optimize.curve_fit(
            _linear_model,
            x_fit,
            y_fit,
            p0=p0_lin,
            bounds=bounds_lin,
            maxfev=10000,
            jac=_linear_model_jac,
        )

        a_lin, b_lin = lin_popt
        lin_alpha = float(a_lin + b_lin)

        ci_low_lin, ci_high_lin = _compute_ci_monte_carlo(
            lin_popt, lin_pcov, _linear_model, ci_replicates
        )
        lin_success = True

    except (RuntimeError, ValueError) as exc:
        logger.debug("Linear model fit failed: %s", exc)

    # Select best model
    # Following asymptoticMK R package: prefer exponential unless CI > 100, then use linear
    # Also compare AIC when both succeed
    use_exponential = False

    if exp_success and lin_success:
        # Calculate AIC for both models
        exp_residuals = y_fit - _exponential_model(x_fit, *exp_popt)
        lin_residuals = y_fit - _linear_model(x_fit, *lin_popt)
        exp_rss = float(np.sum(exp_residuals**2))
        lin_rss = float(np.sum(lin_residuals**2))

        exp_aic = _compute_aic(len(x_fit), exp_rss, 3)
        lin_aic = _compute_aic(len(x_fit), lin_rss, 2)

        # Following asymptoticMK: use linear if exponential CI exceeds 100
        # or if linear AIC is better
        if exp_ci_width > 100:
            use_exponential = False
        elif lin_aic < exp_aic:
            use_exponential = False
        else:
            use_exponential = True
    elif exp_success:
        use_exponential = exp_ci_width <= 100
    elif lin_success:
        use_exponential = False
    else:
        # Both failed, use last value within the fit window (not the
        # untrimmed data -- both attempted fits were against y_fit).
        result = AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_x_values=valid_centers,
            alpha_asymptotic=y_fit[-1],
            ci_low=y_fit[-1],
            ci_high=y_fit[-1],
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            pn_total=pn_total,
            ps_total=ps_total,
            ci_method=ci_method,
            sfs_mode=sfs_mode,
        )
        return _attach_omega(result, agg.ln_total, agg.ls_total)

    if ci_method == "bootstrap":
        # Hold to the model the AIC step selected so the CI is comparable to MC.
        pooled = [p for g in gene_data for p in g.polymorphisms]
        if use_exponential:
            ci_low_exp, ci_high_exp = _compute_ci_bootstrap(
                pooled,
                dn,
                ds,
                bin_edges,
                bin_centers,
                num_bins,
                _exponential_model,
                p0_exp,
                bounds_exp,
                frequency_cutoffs,
                exp_alpha,
                n_replicates=ci_replicates,
                seed=seed,
                sfs_mode=sfs_mode,
                workers=workers,
                model_jac=_exponential_model_jac,
            )
        else:
            ci_low_lin, ci_high_lin = _compute_ci_bootstrap(
                pooled,
                dn,
                ds,
                bin_edges,
                bin_centers,
                num_bins,
                _linear_model,
                p0_lin,
                bounds_lin,
                frequency_cutoffs,
                lin_alpha,
                n_replicates=ci_replicates,
                seed=seed,
                sfs_mode=sfs_mode,
                workers=workers,
                model_jac=_linear_model_jac,
            )

    if use_exponential:
        a, b, c = exp_popt
        result = AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_x_values=valid_centers,
            alpha_asymptotic=exp_alpha,
            ci_low=ci_low_exp,
            ci_high=ci_high_exp,
            fit_a=float(a),
            fit_b=float(b),
            fit_c=float(c),
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            model_type="exponential",
            pn_total=pn_total,
            ps_total=ps_total,
            ci_method=ci_method,
            sfs_mode=sfs_mode,
        )
    else:
        a_lin, b_lin = lin_popt
        result = AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_x_values=valid_centers,
            alpha_asymptotic=lin_alpha,
            ci_low=ci_low_lin,
            ci_high=ci_high_lin,
            fit_a=float(a_lin),
            fit_b=float(b_lin),
            fit_c=0.0,
            dn=dn,
            ds=ds,
            num_genes=agg.num_genes,
            model_type="linear",
            pn_total=pn_total,
            ps_total=ps_total,
            ci_method=ci_method,
            sfs_mode=sfs_mode,
        )
    return _attach_omega(result, agg.ln_total, agg.ls_total)


def asymptotic_mk_test(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    num_bins: int = 10,
    bootstrap_replicates: int = 100,
    frequency_cutoffs: tuple[float, float] = (0.1, 0.9),
    pool_polymorphisms: bool = False,
    sfs_mode: str = "at",
    workers: int = 1,
) -> AsymptoticMKResult:
    """Perform the asymptotic McDonald-Kreitman test.

    This method accounts for slightly deleterious mutations by examining
    how alpha changes across the frequency spectrum and extrapolating
    to the asymptotic value at x=1.

    Based on Messer & Petrov (2013) PNAS.

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup: SequenceSet or path to FASTA file for outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        num_bins: Number of frequency bins (default 10)
        bootstrap_replicates: Number of bootstrap replicates for CI (default 100)
        frequency_cutoffs: (low, high) frequency range for fitting (default 0.1-0.9)
        pool_polymorphisms: If True, consider sites polymorphic in either
            population (libsequence convention). Frequencies are still
            calculated from ingroup only. If False (default), only consider
            ingroup polymorphisms (DnaSP/original MK convention).

    Returns:
        AsymptoticMKResult containing test statistics
    """
    _validate_sfs_mode(sfs_mode)

    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Create aligned pair
    pair = AlignedPair(ingroup=ingroup, outgroup=outgroup, genetic_code=code)
    ln_total, ls_total = pair.count_total_sites()

    # Count divergence
    dn = 0
    ds = 0

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            result = pair.classify_fixed_difference(codon_idx)
            if result is not None:
                nonsyn, syn = result
                dn += nonsyn
                ds += syn

    # Collect polymorphisms with their frequencies
    poly_data: list[tuple[float, str]] = []  # (frequency, type: 'N' or 'S')

    # Get polymorphic sites based on mode
    if pool_polymorphisms:
        poly_sites = pair.polymorphic_sites_pooled()
    else:
        poly_sites = pair.polymorphic_sites_ingroup()

    for codon_idx in poly_sites:
        codons = list(ingroup.codon_set_clean(codon_idx))
        if len(codons) < 2:
            continue

        # Polarizing needs an allele the ingroup shares with the outgroup.
        state = ingroup.derived_state(outgroup, codon_idx)
        if state is None:
            continue
        derived_freq, _ = state

        if derived_freq <= 0 or derived_freq >= 1:
            continue

        # Classify the polymorphism
        result = pair.classify_polymorphism(codon_idx)
        if result is not None:
            nonsyn, syn = result
            for _ in range(nonsyn):
                poly_data.append((derived_freq, "N"))
            for _ in range(syn):
                poly_data.append((derived_freq, "S"))

    # Create frequency bins
    bin_edges = _frequency_bin_edges(num_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Bin polymorphisms by derived allele frequency
    pn_per_bin = np.zeros(num_bins)
    ps_per_bin = np.zeros(num_bins)

    for freq, poly_type in poly_data:
        # Find which bin this frequency falls into
        bin_idx = int(np.searchsorted(bin_edges[1:], freq, side="right"))
        bin_idx = min(bin_idx, num_bins - 1)
        if poly_type == "N":
            pn_per_bin[bin_idx] += 1
        else:
            ps_per_bin[bin_idx] += 1

    pn_per_bin, ps_per_bin = _apply_sfs_mode(pn_per_bin, ps_per_bin, sfs_mode)

    # Calculate per-bin alpha (matching asymptoticMK R package under "at" mode;
    # under "above" mode this is α(>x) on the cumulative SFS).
    alpha_values = []
    valid_bins = []
    valid_centers = []

    for i in range(num_bins):
        pn_bin = float(pn_per_bin[i])
        ps_bin = float(ps_per_bin[i])

        if ds > 0 and dn > 0 and ps_bin > 0:
            alpha_x = 1.0 - (ds / dn) * (pn_bin / ps_bin)
            alpha_values.append(alpha_x)
            valid_bins.append(i)
            valid_centers.append(bin_centers[i])

    if len(valid_centers) < 3:
        # Not enough data points for curve fitting
        simple_alpha = None
        if ds > 0 and dn > 0:
            total_pn = sum(1 for _, t in poly_data if t == "N")
            total_ps = sum(1 for _, t in poly_data if t == "S")
            if total_ps > 0:
                simple_alpha = 1.0 - (ds * total_pn) / (dn * total_ps)

        result = AsymptoticMKResult(
            frequency_bins=list(bin_centers),
            alpha_by_freq=alpha_values,
            alpha_x_values=valid_centers,
            alpha_asymptotic=simple_alpha or 0.0,
            ci_low=simple_alpha or 0.0,
            ci_high=simple_alpha or 0.0,
            dn=dn,
            ds=ds,
            ci_method="bootstrap",
            sfs_mode=sfs_mode,
        )
        return _attach_omega(result, ln_total, ls_total)

    # Fit exponential model: α(x) = a + b * exp(-c * x)
    x_data = np.array(valid_centers)
    y_data = np.array(alpha_values)
    x_fit, y_fit = _mask_to_frequency_cutoffs(x_data, y_data, frequency_cutoffs)

    try:
        # Initial guesses: a is asymptotic value, b is negative for increasing alpha
        p0 = [y_fit[-1], y_fit[0] - y_fit[-1], 5.0]

        # Bounds: b can be negative (typical for increasing alpha with frequency)
        bounds = ([-2.0, -2.0, 0.001], [2.0, 2.0, 100.0])

        popt, _ = optimize.curve_fit(
            _exponential_model,
            x_fit,
            y_fit,
            p0=p0,
            bounds=bounds,
            maxfev=10000,
            jac=_exponential_model_jac,
        )

        a, b, c = popt
        # Evaluate the model at x=1: α(1) = a + b * exp(-c)
        alpha_asymp = float(a + b * np.exp(-c))

    except (RuntimeError, ValueError):
        # Curve fitting failed, use last value within the fit window (not the
        # untrimmed data -- the attempted fit itself was against y_fit).
        a, b, c = y_fit[-1], 0.0, 1.0
        alpha_asymp = y_fit[-1]

    # Shares the case-resampling bootstrap with the aggregated path.
    # Replicates whose curve_fit fails are dropped (rather than imputed
    # with the last bin's alpha as the legacy per-gene path did); CI
    # degenerates to the point estimate when fewer than half succeed.
    ci_low, ci_high = _compute_ci_bootstrap(
        poly_data,
        dn,
        ds,
        bin_edges,
        bin_centers,
        num_bins,
        _exponential_model,
        [float(a), float(b), float(c)],
        bounds,
        frequency_cutoffs,
        alpha_asymp,
        n_replicates=bootstrap_replicates,
        seed=42,
        sfs_mode=sfs_mode,
        workers=workers,
        model_jac=_exponential_model_jac,
    )

    result = AsymptoticMKResult(
        frequency_bins=list(bin_centers),
        alpha_by_freq=alpha_values,
        alpha_x_values=valid_centers,
        alpha_asymptotic=alpha_asymp,
        ci_low=ci_low,
        ci_high=ci_high,
        fit_a=float(a),
        fit_b=float(b),
        fit_c=float(c),
        dn=dn,
        ds=ds,
        ci_method="bootstrap",
        sfs_mode=sfs_mode,
    )
    return _attach_omega(result, ln_total, ls_total)
