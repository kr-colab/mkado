"""Tarone-Greenland alpha (α_TG) implementation.

Implements the weighted estimator from:
    Stoletzki N, Eyre-Walker A (2011) Estimation of the Neutrality Index.
    Mol Biol Evol 28(1):63-70. doi:10.1093/molbev/msq249

The NI_TG formula (Equation 3):

         Σᵢ (Dₛᵢ × Pₙᵢ) / (Pₛᵢ + Dₛᵢ)
NI_TG = ────────────────────────────────
         Σᵢ (Dₙᵢ × Pₛᵢ) / (Pₛᵢ + Dₛᵢ)

And α_TG = 1 - NI_TG (proportion of adaptive nonsynonymous substitutions).

The weighting by 1/(Pₛᵢ + Dₛᵢ) provides an unbiased estimate even with
heterogeneity across genes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mkado.analysis.asymptotic import PolymorphismData, sum_site_totals
from mkado.analysis.statistics import omega_decomposition


@dataclass
class AlphaTGResult:
    """Results from Tarone-Greenland alpha estimation.

    Per-gene ``(Dn, Ds, Pn, Ps)`` counts follow the standard MK definitions
    (see :class:`mkado.analysis.mk_test.MKResult`); ``min_frequency``
    filtering applied to each gene's polymorphism extraction propagates
    into the totals reported here.
    """

    alpha_tg: float
    """Proportion of adaptive substitutions (1 - NI_TG)."""

    ni_tg: float
    """The underlying weighted neutrality index."""

    ci_low: float | None
    """Lower bound of 95% bootstrap confidence interval for alpha_tg.

    ``None`` when fewer than two genes are available: the gene-resampling
    bootstrap always redraws the same single gene, so no interval can be
    estimated.
    """

    ci_high: float | None
    """Upper bound of 95% bootstrap confidence interval for alpha_tg.

    ``None`` under the same one-gene condition as ``ci_low``.
    """

    num_genes: int
    """Number of genes used in the calculation."""

    dn_total: int
    """Total nonsynonymous divergence across all genes."""

    ds_total: int
    """Total synonymous divergence across all genes."""

    pn_total: int
    """Total nonsynonymous polymorphism across all genes."""

    ps_total: int
    """Total synonymous polymorphism across all genes."""

    ln: float | None = None
    """Total non-synonymous sites (Nei-Gojobori) summed across genes."""
    ls: float | None = None
    """Total synonymous sites (Nei-Gojobori) summed across genes."""
    omega: float | None = None
    """``(Dn_total/Ds_total) * (Ls/Ln)`` — pooled dN/dS."""
    omega_a: float | None = None
    """Adaptive rate ``alpha_tg * omega`` (Gossmann, Keightley & Eyre-Walker 2012)."""
    omega_na: float | None = None
    """Non-adaptive component ``(1 - alpha_tg) * omega``."""
    omega_ci_low: float | None = None
    """Lower 95% bootstrap CI for ``omega`` (gene-level resampling varies Dn, Ds, Ln, Ls)."""
    omega_ci_high: float | None = None
    """Upper 95% bootstrap CI for ``omega``."""
    omega_a_ci_low: float | None = None
    """Lower 95% bootstrap CI for ``omega_a``."""
    omega_a_ci_high: float | None = None
    """Upper 95% bootstrap CI for ``omega_a``."""
    omega_na_ci_low: float | None = None
    """Lower 95% bootstrap CI for ``omega_na``."""
    omega_na_ci_high: float | None = None
    """Upper 95% bootstrap CI for ``omega_na``."""
    ci_method: str = "bootstrap"
    """Method used for the CI. Always ``"bootstrap"`` (gene-level resampling)
    for ``AlphaTGResult`` — the weighted estimator has no parametric MC analog."""

    def __str__(self) -> str:
        """Return a human-readable string representation."""
        ci_low_str = f"{self.ci_low:.4f}" if self.ci_low is not None else "NA"
        ci_high_str = f"{self.ci_high:.4f}" if self.ci_high is not None else "NA"
        lines = [
            "Tarone-Greenland Alpha (Stoletzki & Eyre-Walker 2011):",
            f"  α_TG:  {self.alpha_tg:.4f} "
            f"(95% CI [{self.ci_method}]: {ci_low_str} - {ci_high_str})",
            f"  NI_TG: {self.ni_tg:.4f}",
            f"  Divergence:    Dn={self.dn_total}, Ds={self.ds_total}",
            f"  Polymorphism:  Pn={self.pn_total}, Ps={self.ps_total}",
            f"  Genes: {self.num_genes}",
        ]
        if self.ln is not None and self.ls is not None:
            lines.append(f"  Sites:         Ln={self.ln:.2f}, Ls={self.ls:.2f}")
        if self.omega is not None:
            omega_a_str = f"{self.omega_a:.4f}" if self.omega_a is not None else "NA"
            omega_na_str = f"{self.omega_na:.4f}" if self.omega_na is not None else "NA"
            lines.append(
                f"  omega:         {self.omega:.4f} "
                f"(omega_a={omega_a_str}, omega_na={omega_na_str})"
            )
            if self.omega_ci_low is not None and self.omega_ci_high is not None:
                lines.append(
                    f"    omega 95% CI:    ({self.omega_ci_low:.4f}, {self.omega_ci_high:.4f})"
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
        return "\n".join(lines)

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "alpha_tg": self.alpha_tg,
            "ni_tg": self.ni_tg,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "num_genes": self.num_genes,
            "dn_total": self.dn_total,
            "ds_total": self.ds_total,
            "pn_total": self.pn_total,
            "ps_total": self.ps_total,
            "ln": self.ln,
            "ls": self.ls,
            "omega": self.omega,
            "omega_a": self.omega_a,
            "omega_na": self.omega_na,
            "omega_ci_low": self.omega_ci_low,
            "omega_ci_high": self.omega_ci_high,
            "omega_a_ci_low": self.omega_a_ci_low,
            "omega_a_ci_high": self.omega_a_ci_high,
            "omega_na_ci_low": self.omega_na_ci_low,
            "omega_na_ci_high": self.omega_na_ci_high,
            "ci_method": self.ci_method,
        }


def compute_ni_tg(gene_data: list[PolymorphismData]) -> float | None:
    """Compute NI_TG (Equation 3 from Stoletzki & Eyre-Walker 2011).

    NI_TG = Σᵢ (Dₛᵢ × Pₙᵢ) / (Pₛᵢ + Dₛᵢ)
            ─────────────────────────────
            Σᵢ (Dₙᵢ × Pₛᵢ) / (Pₛᵢ + Dₛᵢ)

    Args:
        gene_data: List of PolymorphismData from individual genes.

    Returns:
        The weighted neutrality index, or None if denominator is zero.
    """
    numerator = 0.0
    denominator = 0.0

    for gene in gene_data:
        # Count Pn and Ps for this gene
        pn = sum(1 for _, ptype in gene.polymorphisms if ptype == "N")
        ps = sum(1 for _, ptype in gene.polymorphisms if ptype == "S")

        # Weight factor: 1 / (Ps + Ds)
        weight = ps + gene.ds
        if weight > 0:
            numerator += (gene.ds * pn) / weight
            denominator += (gene.dn * ps) / weight

    if denominator <= 0:
        return None

    return numerator / denominator


def alpha_tg_from_gene_data(
    gene_data: list[PolymorphismData],
    bootstrap_replicates: int = 1000,
    seed: int | None = None,
) -> AlphaTGResult:
    """Compute α_TG with bootstrap confidence intervals.

    Args:
        gene_data: List of PolymorphismData from individual genes.
        bootstrap_replicates: Number of bootstrap replicates for CI estimation.
        seed: Random seed for reproducibility (optional).

    Returns:
        AlphaTGResult containing α_TG, NI_TG, and 95% bootstrap CI.
    """
    # Calculate totals
    dn_total = sum(g.dn for g in gene_data)
    ds_total = sum(g.ds for g in gene_data)
    pn_total = sum(sum(1 for _, ptype in g.polymorphisms if ptype == "N") for g in gene_data)
    ps_total = sum(sum(1 for _, ptype in g.polymorphisms if ptype == "S") for g in gene_data)
    num_genes = len(gene_data)

    ln_total, ls_total = sum_site_totals(gene_data)

    # Compute point estimate
    ni_tg = compute_ni_tg(gene_data)

    if ni_tg is None:
        # Cannot compute - return zero with NA-like CI
        return AlphaTGResult(
            alpha_tg=0.0,
            ni_tg=0.0,
            ci_low=0.0,
            ci_high=0.0,
            num_genes=num_genes,
            dn_total=dn_total,
            ds_total=ds_total,
            pn_total=pn_total,
            ps_total=ps_total,
            ln=ln_total,
            ls=ls_total,
        )

    alpha_tg = 1.0 - ni_tg

    # Bootstrap for confidence intervals; resamples genes, so Dn/Ds/Ln/Ls
    # totals also vary per replicate, giving omega its own sampling distribution.
    # With one gene every replicate redraws that same gene, so the interval
    # cannot be estimated; skip the resampling and report it as undefined
    # rather than a zero-width interval that looks like real precision.
    bootstrap_alphas: list[float] = []
    # omega_decomposition returns all three or none, so these stay in lockstep.
    bootstrap_omegas: list[float] = []
    bootstrap_omega_a: list[float] = []
    bootstrap_omega_na: list[float] = []

    if num_genes >= 2:
        rng = np.random.default_rng(seed)
        for _ in range(bootstrap_replicates):
            indices = rng.integers(0, num_genes, size=num_genes)
            boot_data = [gene_data[i] for i in indices]

            boot_ni = compute_ni_tg(boot_data)
            if boot_ni is None:
                continue
            boot_alpha = 1.0 - boot_ni
            bootstrap_alphas.append(boot_alpha)

            boot_dn = sum(g.dn for g in boot_data)
            boot_ds = sum(g.ds for g in boot_data)
            boot_ln, boot_ls = sum_site_totals(boot_data)
            boot_omega, boot_oa, boot_ona = omega_decomposition(
                boot_dn, boot_ds, boot_ln, boot_ls, boot_alpha
            )
            if boot_omega is not None:
                bootstrap_omegas.append(boot_omega)
                bootstrap_omega_a.append(boot_oa)
                bootstrap_omega_na.append(boot_ona)

    if bootstrap_alphas:
        ci_low, ci_high = np.percentile(bootstrap_alphas, [2.5, 97.5])
        ci_low, ci_high = float(ci_low), float(ci_high)
    else:
        # Undefined, not a fabricated zero-width interval: either num_genes < 2
        # (skipped above) or every replicate failed (denominator always zero).
        ci_low = ci_high = None

    if bootstrap_omegas:
        cis = np.percentile(
            [bootstrap_omegas, bootstrap_omega_a, bootstrap_omega_na],
            [2.5, 97.5],
            axis=1,
        )
        omega_ci_low, omega_a_ci_low, omega_na_ci_low = (float(x) for x in cis[0])
        omega_ci_high, omega_a_ci_high, omega_na_ci_high = (float(x) for x in cis[1])
    else:
        omega_ci_low = omega_ci_high = None
        omega_a_ci_low = omega_a_ci_high = None
        omega_na_ci_low = omega_na_ci_high = None

    omega, omega_a, omega_na = omega_decomposition(dn_total, ds_total, ln_total, ls_total, alpha_tg)

    return AlphaTGResult(
        alpha_tg=alpha_tg,
        ni_tg=ni_tg,
        ci_low=ci_low,
        ci_high=ci_high,
        num_genes=num_genes,
        dn_total=dn_total,
        ds_total=ds_total,
        pn_total=pn_total,
        ps_total=ps_total,
        ln=ln_total,
        ls=ls_total,
        omega=omega,
        omega_a=omega_a,
        omega_na=omega_na,
        omega_ci_low=omega_ci_low,
        omega_ci_high=omega_ci_high,
        omega_a_ci_low=omega_a_ci_low,
        omega_a_ci_high=omega_a_ci_high,
        omega_na_ci_low=omega_na_ci_low,
        omega_na_ci_high=omega_na_ci_high,
    )
