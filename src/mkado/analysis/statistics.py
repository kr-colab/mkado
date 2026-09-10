"""Statistical functions for MK test analysis."""

from __future__ import annotations

import math

from scipy import stats


def fishers_exact(dn: int, ds: int, pn: int, ps: int, alternative: str = "two-sided") -> float:
    """Perform Fisher's exact test on a 2x2 contingency table.

    The table is:
                  | Non-synonymous | Synonymous
        -----------------------------------------------
        Divergence |      dn       |     ds
        Polymorphism |    pn       |     ps

    Args:
        dn: Non-synonymous divergence (fixed differences)
        ds: Synonymous divergence (fixed differences)
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms
        alternative: 'two-sided', 'less', or 'greater'

    Returns:
        p-value from Fisher's exact test
    """
    table = [[dn, ds], [pn, ps]]
    _, p_value = stats.fisher_exact(table, alternative=alternative)
    return float(p_value)


def neutrality_index(dn: int, ds: int, pn: int, ps: int) -> float | None:
    """Calculate the Neutrality Index (NI).

    NI = (Pn/Ps) / (Dn/Ds)

    Under neutrality, NI = 1.
    NI > 1 suggests segregating weakly deleterious polymorphisms (excess polymorphism).
    NI < 1 suggests positive selection (excess divergence).

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        Neutrality Index, or None if cannot be calculated (division by zero)
    """
    if ds == 0 or ps == 0 or dn == 0:
        return None

    return (pn / ps) / (dn / ds)


def alpha(dn: int, ds: int, pn: int, ps: int) -> float | None:
    """Calculate alpha, the proportion of adaptive substitutions.

    α = 1 - (Ds × Pn) / (Dn × Ps)

    This per-gene formulation is from Smith and Eyre-Walker (2002):
    https://pubmed.ncbi.nlm.nih.gov/11875568/

    Under neutrality, α = 0.
    α > 0 suggests positive selection.
    α < 0 suggests segregating weakly deleterious polymorphisms.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        Alpha value, or None if cannot be calculated
    """
    if dn == 0 or ps == 0:
        return None

    return 1.0 - (ds * pn) / (dn * ps)


def dos(dn: int, ds: int, pn: int, ps: int) -> float | None:
    """Calculate the Direction of Selection (DoS) statistic.

    DoS = Dn/(Dn+Ds) - Pn/(Pn+Ps)

    From Stoletzki & Eyre-Walker (2011).

    Interpretation:
    - DoS = 0: Neutral evolution
    - DoS > 0: Positive selection (excess adaptive substitutions)
    - DoS < 0: Slightly deleterious polymorphisms

    Advantages over NI: Symmetric around 0, bounded [-1, +1], well-behaved
    with small counts.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        DoS value, or None if cannot be calculated (all counts zero)
    """
    dn_total = dn + ds
    pn_total = pn + ps

    if dn_total == 0 and pn_total == 0:
        return None

    dn_ratio = dn / dn_total if dn_total > 0 else 0.0
    pn_ratio = pn / pn_total if pn_total > 0 else 0.0

    return dn_ratio - pn_ratio


def g_test(dn: int, ds: int, pn: int, ps: int) -> float:
    """Perform G-test (log-likelihood ratio test) on MK table.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms

    Returns:
        p-value from G-test
    """
    observed = [[dn, ds], [pn, ps]]
    total = dn + ds + pn + ps
    if total == 0:
        return 1.0

    row_totals = [dn + ds, pn + ps]
    col_totals = [dn + pn, ds + ps]

    g_stat = 0.0
    for i in range(2):
        for j in range(2):
            obs = observed[i][j]
            if obs > 0:
                expected = row_totals[i] * col_totals[j] / total
                if expected > 0:
                    g_stat += 2 * obs * math.log(obs / expected)

    # G-test statistic follows chi-square with df=1. The survival function
    # keeps precision in the far tail, where 1 - cdf rounds to zero.
    p_value = stats.chi2.sf(g_stat, df=1)
    return float(p_value)


def omega_decomposition(
    dn: int,
    ds: int,
    ln: float | None,
    ls: float | None,
    alpha_value: float | None,
) -> tuple[float | None, float | None, float | None]:
    """Compute omega and its adaptive/non-adaptive decomposition.

    Decomposition introduced by Gossmann, Keightley & Eyre-Walker (2012,
    *Genome Biol Evol* 4(5):658-667) and applied to MK-style counts by
    Coronado-Zamora et al. (2019, *Genome Biol Evol* 11(5):1463-1482)::

        omega    = (Dn / Ds) * (Ls / Ln)   (= dN/dS using N-G site counts)
        omega_a  = alpha * omega
        omega_na = (1 - alpha) * omega = omega - omega_a

    Returns ``(None, None, None)`` if site counts are missing or any
    denominator is zero. ``omega_a`` and ``omega_na`` are ``None``
    when ``alpha_value`` is ``None``.
    """
    if ln is None or ls is None or ln <= 0 or ls <= 0 or ds == 0:
        return (None, None, None)

    omega = (dn * ls) / (ds * ln)
    if alpha_value is None:
        return (omega, None, None)
    return (omega, alpha_value * omega, (1.0 - alpha_value) * omega)


def confidence_interval_alpha(
    dn: int, ds: int, pn: int, ps: int, confidence: float = 0.95
) -> tuple[float, float] | None:
    """Calculate confidence interval for alpha using bootstrap-like approximation.

    Uses the delta method to approximate the variance of alpha.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms
        confidence: Confidence level (default 0.95)

    Returns:
        Tuple of (lower_bound, upper_bound), or None if cannot be calculated
    """
    a = alpha(dn, ds, pn, ps)
    if a is None:
        return None

    # Use Wilson score interval approximation for small counts
    # For simplicity, use normal approximation with pseudo-counts
    dn_p = max(dn, 0.5)
    ds_p = max(ds, 0.5)
    pn_p = max(pn, 0.5)
    ps_p = max(ps, 0.5)

    # Variance approximation using delta method
    # alpha = 1 - (ds*pn)/(dn*ps)
    # Let r = (ds*pn)/(dn*ps)
    r = (ds_p * pn_p) / (dn_p * ps_p)

    # Var(r) ≈ r^2 * (1/dn + 1/ds + 1/pn + 1/ps) using delta method
    var_r = r**2 * (1 / dn_p + 1 / ds_p + 1 / pn_p + 1 / ps_p)

    # Var(alpha) = Var(r) since alpha = 1 - r
    se_alpha = math.sqrt(var_r)

    # z-score for confidence level
    z = stats.norm.ppf((1 + confidence) / 2)

    lower = a - z * se_alpha
    upper = a + z * se_alpha

    return (lower, upper)
