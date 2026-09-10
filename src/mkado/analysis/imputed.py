"""Imputed McDonald-Kreitman test implementation.

Based on Murga-Moreno et al. (2022) G3: Genes|Genomes|Genetics.
https://doi.org/10.1093/g3journal/jkac206

The impMKT corrects for slightly deleterious mutations by imputing the number
of weakly deleterious nonsynonymous polymorphisms segregating at low frequency,
rather than discarding all low-frequency variants. This retains more data and
increases power to detect positive selection at the gene level.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from mkado.analysis.asymptotic import PolymorphismData, sum_site_totals
from mkado.analysis.statistics import fishers_exact, omega_decomposition

logger = logging.getLogger(__name__)


@dataclass
class ImputedMKResult:
    """Results from an imputed McDonald-Kreitman test.

    The imputed test partitions ``Pn`` into a weakly-deleterious fraction
    (``Pwd``) — imputed from the synonymous SFS scaled by the observed
    Pn/Ps ratio at low derived allele frequencies — and a neutral
    fraction (``Pn_neutral = Pn_total - Pwd``). ``alpha`` is computed
    against ``Pn_neutral`` rather than the raw ``Pn`` count, removing the
    downward bias caused by slightly deleterious nonsynonymous variants.

    ``Dn`` and ``Ds`` follow the standard MK definitions. ``Pn_total`` and
    ``Ps_total`` are the unfiltered ingroup polymorphism counts (the imputed
    test uses the entire SFS, including singletons).
    """

    alpha: float | None
    """Proportion of adaptive substitutions, computed against ``Pn_neutral``."""
    p_value: float
    """Fisher's exact test p-value on the imputed (Dn, Ds, Pn_neutral, Ps_total) table."""
    pn_neutral: float
    """Estimated neutral nonsynonymous polymorphism count (``Pn_total - Pwd``)."""
    pwd: float
    """Imputed count of weakly-deleterious nonsynonymous polymorphisms."""
    dn: int
    """Non-synonymous fixed differences."""
    ds: int
    """Synonymous fixed differences."""
    pn_total: int
    """Total non-synonymous ingroup polymorphisms (no frequency filter applied)."""
    ps_total: int
    """Total synonymous ingroup polymorphisms (no frequency filter applied)."""
    d: float | None = None
    """Estimated DFE deleterious fraction (Murga-Moreno et al. 2022, eq. 3)."""
    b: float | None = None
    """Estimated DFE weakly-deleterious fraction."""
    f: float | None = None
    """Estimated DFE neutral fraction."""
    cutoff: float = 0.15
    """Derived allele frequency below which polymorphisms are treated as the imputation pool."""
    ln: float | None = None
    """Total non-synonymous sites (Nei-Gojobori) over analyzed codons."""
    ls: float | None = None
    """Total synonymous sites (Nei-Gojobori) over analyzed codons."""
    omega: float | None = None
    """``(Dn/Ds) * (Ls/Ln)`` — dN/dS ratio."""
    omega_a: float | None = None
    """Adaptive rate ``alpha * omega`` with imputed alpha (Gossmann, Keightley & Eyre-Walker 2012)."""
    omega_na: float | None = None
    """Non-adaptive component ``(1 - alpha) * omega``."""
    alpha_ci_low: float | None = None
    """Lower 95% bootstrap CI for ``alpha``. ``None`` when ``n_bootstrap=0``."""
    alpha_ci_high: float | None = None
    """Upper 95% bootstrap CI for ``alpha``."""
    ci_method: str | None = None
    """``"bootstrap"`` when CI was computed, ``None`` when ``n_bootstrap=0``."""
    stop_blocked_pairs: int = 0
    """Codon pairs with no stop-free mutational ordering, dropped from Dn/Ds/Pn/Ps."""

    # omega_a/omega_na CIs are pure arithmetic transforms of the alpha CI scaled
    # by the (constant) point-estimate omega; mirrors AsymptoticMKResult.
    @property
    def omega_a_ci_low(self) -> float | None:
        """Lower 95% CI for ``omega_a`` = ``alpha_ci_low * omega``."""
        if self.omega is None or self.alpha_ci_low is None:
            return None
        return self.alpha_ci_low * self.omega

    @property
    def omega_a_ci_high(self) -> float | None:
        """Upper 95% CI for ``omega_a`` = ``alpha_ci_high * omega``."""
        if self.omega is None or self.alpha_ci_high is None:
            return None
        return self.alpha_ci_high * self.omega

    @property
    def omega_na_ci_low(self) -> float | None:
        """Lower 95% CI for ``omega_na`` = ``(1 - alpha_ci_high) * omega`` (flipped)."""
        if self.omega is None or self.alpha_ci_high is None:
            return None
        return (1.0 - self.alpha_ci_high) * self.omega

    @property
    def omega_na_ci_high(self) -> float | None:
        """Upper 95% CI for ``omega_na`` = ``(1 - alpha_ci_low) * omega``."""
        if self.omega is None or self.alpha_ci_low is None:
            return None
        return (1.0 - self.alpha_ci_low) * self.omega

    def __str__(self) -> str:
        alpha_str = f"{self.alpha:.4f}" if self.alpha is not None else "NA"
        lines = [
            "Imputed MK Test Results:",
            f"  Divergence:    Dn={self.dn}, Ds={self.ds}",
            f"  Polymorphism:  Pn={self.pn_total}, Ps={self.ps_total}",
            f"  DAF cutoff:    {self.cutoff}",
            f"  Imputed Pwd:   {self.pwd:.2f}",
            f"  Pn (neutral):  {self.pn_neutral:.2f}",
            f"  Alpha (α):     {alpha_str}",
            f"  p-value:       {self.p_value:.4g}",
        ]
        if self.alpha_ci_low is not None and self.alpha_ci_high is not None:
            lines.append(
                f"  Alpha 95% CI [{self.ci_method}]: "
                f"({self.alpha_ci_low:.4f}, {self.alpha_ci_high:.4f})"
            )
        if self.ln is not None and self.ls is not None:
            lines.append(f"  Sites:         Ln={self.ln:.2f}, Ls={self.ls:.2f}")
        if self.omega is not None:
            omega_a_str = f"{self.omega_a:.4f}" if self.omega_a is not None else "NA"
            omega_na_str = f"{self.omega_na:.4f}" if self.omega_na is not None else "NA"
            lines.append(
                f"  omega:         {self.omega:.4f} "
                f"(omega_a={omega_a_str}, omega_na={omega_na_str})"
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
        if self.d is not None:
            lines.append(f"  DFE fractions: d={self.d:.4f}, b={self.b:.4f}, f={self.f:.4f}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        result = {
            "alpha": self.alpha,
            "p_value": self.p_value,
            "pn_neutral": self.pn_neutral,
            "pwd": self.pwd,
            "dn": self.dn,
            "ds": self.ds,
            "pn_total": self.pn_total,
            "ps_total": self.ps_total,
            "cutoff": self.cutoff,
            "ln": self.ln,
            "ls": self.ls,
            "omega": self.omega,
            "omega_a": self.omega_a,
            "omega_na": self.omega_na,
            "alpha_ci_low": self.alpha_ci_low,
            "alpha_ci_high": self.alpha_ci_high,
            "ci_method": self.ci_method,
            "omega_a_ci_low": self.omega_a_ci_low,
            "omega_a_ci_high": self.omega_a_ci_high,
            "omega_na_ci_low": self.omega_na_ci_low,
            "omega_na_ci_high": self.omega_na_ci_high,
        }
        if self.d is not None:
            result["d"] = self.d
            result["b"] = self.b
            result["f"] = self.f
        return result


def _compute_imputed(
    polymorphisms: list[tuple[float, str]],
    dn: int,
    ds: int,
    cutoff: float,
    num_synonymous_sites: float | None,
    num_nonsynonymous_sites: float | None,
) -> ImputedMKResult:
    """Core imputed MK algorithm on pooled data."""
    # Split polymorphisms at cutoff
    pn_low = sum(1 for freq, t in polymorphisms if t == "N" and freq <= cutoff)
    pn_high = sum(1 for freq, t in polymorphisms if t == "N" and freq > cutoff)
    ps_low = sum(1 for freq, t in polymorphisms if t == "S" and freq <= cutoff)
    ps_high = sum(1 for freq, t in polymorphisms if t == "S" and freq > cutoff)

    pn_total = pn_low + pn_high
    ps_total = ps_low + ps_high

    # Compute ratio of low to high frequency synonymous polymorphisms
    if ps_high == 0:
        ratio_ps = 0.0
    else:
        ratio_ps = ps_low / ps_high

    # Impute weakly deleterious nonsynonymous polymorphisms
    pwd = pn_low - pn_high * ratio_ps
    pwd = max(pwd, 0.0)

    # Neutral nonsynonymous polymorphisms
    pn_neutral = pn_total - pwd

    # Corrected alpha
    if dn == 0 or ps_total == 0:
        alpha_val = None
    else:
        alpha_val = 1.0 - (pn_neutral / ps_total) * (ds / dn)

    # Fisher's exact test on corrected table
    pn_neutral_int = round(pn_neutral)
    p_value = fishers_exact(dn, ds, pn_neutral_int, ps_total)

    # DFE fractions
    d_frac = None
    b_frac = None
    f_frac = None

    if num_synonymous_sites is not None and num_nonsynonymous_sites is not None:
        m0 = num_synonymous_sites
        mi = num_nonsynonymous_sites
        if mi > 0 and ps_total > 0:
            b_frac = (pwd / ps_total) * (m0 / mi)
            f_frac = (m0 * pn_neutral) / (mi * ps_total)
            d_frac = 1.0 - f_frac - b_frac

    # ``num_*_sites`` are aliases for Ls (m0) and Ln (mi) — the DFE inputs
    # double as Nei-Gojobori site totals for the omega decomposition.
    omega, omega_a, omega_na = omega_decomposition(
        dn, ds, num_nonsynonymous_sites, num_synonymous_sites, alpha_val
    )

    return ImputedMKResult(
        alpha=alpha_val,
        p_value=p_value,
        pn_neutral=pn_neutral,
        pwd=pwd,
        dn=dn,
        ds=ds,
        pn_total=pn_total,
        ps_total=ps_total,
        d=d_frac,
        b=b_frac,
        f=f_frac,
        cutoff=cutoff,
        ln=num_nonsynonymous_sites,
        ls=num_synonymous_sites,
        omega=omega,
        omega_a=omega_a,
        omega_na=omega_na,
    )


def _bootstrap_imputed_alpha(
    polymorphisms: list[tuple[float, str]],
    dn: int,
    ds: int,
    cutoff: float,
    num_synonymous_sites: float | None,
    num_nonsynonymous_sites: float | None,
    n_replicates: int,
    seed: int,
) -> tuple[float, float] | None:
    """Bootstrap CI on imputed alpha by case-resampling the polymorphism list.

    Vectorized: pre-extract ``freqs`` / ``is_n`` / ``is_low`` numpy arrays
    once, then per replicate index with ``rng.integers`` and count via
    ``np.count_nonzero``. The CI only needs alpha, so we inline the imputed
    formula and skip Fisher's-exact, DFE fractions, and omega — all of which
    ``_compute_imputed`` recomputes on every call but none of which feed the
    percentile.

    Returns ``None`` if fewer than half the replicates yield a defined alpha.
    """
    n = len(polymorphisms)
    if n == 0 or n_replicates <= 0:
        return None
    if dn == 0:
        # alpha is undefined for every replicate; skip the loop entirely.
        return None

    freqs = np.fromiter((f for f, _ in polymorphisms), dtype=float, count=n)
    is_n = np.fromiter((t == "N" for _, t in polymorphisms), dtype=bool, count=n)
    is_low = freqs <= cutoff

    rng = np.random.default_rng(seed)
    bootstrap_alphas: list[float] = []
    for _ in range(n_replicates):
        sample = rng.integers(0, n, size=n)
        boot_n = is_n[sample]
        boot_low = is_low[sample]
        pn_low = int(np.count_nonzero(boot_n & boot_low))
        pn_high = int(np.count_nonzero(boot_n & ~boot_low))
        ps_low = int(np.count_nonzero(~boot_n & boot_low))
        ps_high = int(np.count_nonzero(~boot_n & ~boot_low))

        ps_total = ps_low + ps_high
        if ps_total == 0:
            continue
        ratio_ps = ps_low / ps_high if ps_high else 0.0
        pwd = max(pn_low - pn_high * ratio_ps, 0.0)
        pn_neutral = (pn_low + pn_high) - pwd
        bootstrap_alphas.append(1.0 - (pn_neutral / ps_total) * (ds / dn))

    if len(bootstrap_alphas) < n_replicates // 2:
        logger.debug(
            "Imputed bootstrap CI: only %d/%d replicates yielded defined alpha; skipping CI",
            len(bootstrap_alphas),
            n_replicates,
        )
        return None

    return (
        float(np.percentile(bootstrap_alphas, 2.5)),
        float(np.percentile(bootstrap_alphas, 97.5)),
    )


def imputed_mk_test(
    gene_data: PolymorphismData,
    cutoff: float = 0.15,
    num_synonymous_sites: float | None = None,
    num_nonsynonymous_sites: float | None = None,
    n_bootstrap: int = 0,
    seed: int = 42,
) -> ImputedMKResult:
    """Perform the imputed McDonald-Kreitman test on a single gene.

    The impMKT corrects for slightly deleterious mutations by imputing the
    count of weakly deleterious nonsynonymous polymorphisms from the ratio
    of low- to high-frequency synonymous polymorphisms.

    Based on Murga-Moreno et al. (2022) G3.
    https://doi.org/10.1093/g3journal/jkac206

    Args:
        gene_data: Polymorphism and divergence data from extract_polymorphism_data()
        cutoff: DAF cutoff separating low and high frequency (default 0.15)
        num_synonymous_sites: Number of synonymous sites (m0) for DFE fractions
        num_nonsynonymous_sites: Number of nonsynonymous sites (mi) for DFE fractions
        n_bootstrap: Number of bootstrap replicates for the alpha CI. ``0`` (default)
            disables CI computation and preserves byte-identical legacy output.
        seed: Random seed for the bootstrap (default 42).

    Returns:
        ImputedMKResult with corrected alpha and optionally DFE fractions
    """
    ls = num_synonymous_sites if num_synonymous_sites is not None else gene_data.ls
    ln = num_nonsynonymous_sites if num_nonsynonymous_sites is not None else gene_data.ln
    result = _compute_imputed(
        gene_data.polymorphisms,
        gene_data.dn,
        gene_data.ds,
        cutoff,
        ls,
        ln,
    )
    result.stop_blocked_pairs = gene_data.stop_blocked_pairs
    if n_bootstrap > 0:
        ci = _bootstrap_imputed_alpha(
            gene_data.polymorphisms,
            gene_data.dn,
            gene_data.ds,
            cutoff,
            ls,
            ln,
            n_bootstrap,
            seed,
        )
        if ci is not None:
            result.alpha_ci_low, result.alpha_ci_high = ci
            result.ci_method = "bootstrap"
    return result


def imputed_mk_test_multi(
    gene_data: list[PolymorphismData],
    cutoff: float = 0.15,
    num_synonymous_sites: float | None = None,
    num_nonsynonymous_sites: float | None = None,
    n_bootstrap: int = 0,
    seed: int = 42,
) -> ImputedMKResult:
    """Perform the imputed MK test on pooled data from multiple genes.

    Pools all polymorphisms and divergence counts across genes, then
    runs the imputed MK algorithm on the combined data.

    Args:
        gene_data: List of PolymorphismData from individual genes
        cutoff: DAF cutoff separating low and high frequency (default 0.15)
        num_synonymous_sites: Number of synonymous sites (m0) for DFE fractions
        num_nonsynonymous_sites: Number of nonsynonymous sites (mi) for DFE fractions
        n_bootstrap: Number of bootstrap replicates for the alpha CI. ``0`` (default)
            disables CI computation and preserves byte-identical legacy output.
        seed: Random seed for the bootstrap (default 42).

    Returns:
        ImputedMKResult with corrected alpha from pooled data
    """
    all_polymorphisms: list[tuple[float, str]] = []
    dn_total = 0
    ds_total = 0
    stop_blocked_total = 0

    for g in gene_data:
        all_polymorphisms.extend(g.polymorphisms)
        dn_total += g.dn
        ds_total += g.ds
        stop_blocked_total += g.stop_blocked_pairs

    ln_agg, ls_agg = sum_site_totals(gene_data)
    ls = num_synonymous_sites if num_synonymous_sites is not None else ls_agg
    ln = num_nonsynonymous_sites if num_nonsynonymous_sites is not None else ln_agg

    result = _compute_imputed(
        all_polymorphisms,
        dn_total,
        ds_total,
        cutoff,
        ls,
        ln,
    )
    result.stop_blocked_pairs = stop_blocked_total
    if n_bootstrap > 0:
        ci = _bootstrap_imputed_alpha(
            all_polymorphisms,
            dn_total,
            ds_total,
            cutoff,
            ls,
            ln,
            n_bootstrap,
            seed,
        )
        if ci is not None:
            result.alpha_ci_low, result.alpha_ci_high = ci
            result.ci_method = "bootstrap"
    return result
