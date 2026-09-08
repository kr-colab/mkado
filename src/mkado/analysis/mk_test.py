"""Standard McDonald-Kreitman test implementation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from mkado.analysis.statistics import (
    alpha,
    dos,
    fishers_exact,
    neutrality_index,
    omega_decomposition,
)
from mkado.core.alignment import AlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class MKResult:
    """Results from a McDonald-Kreitman test.

    Site counts follow the standard MK definitions:

    - **Dn / Ds**: number of fixed differences between ingroup and outgroup
      that change (Dn) or do not change (Ds) the encoded amino acid under
      the supplied genetic code.
    - **Pn / Ps**: number of ingroup polymorphic sites whose derived allele
      causes a non-synonymous (Pn) or synonymous (Ps) substitution. Sites
      are filtered by ``min_frequency`` (default ``0.0``); the
      ``--no-singletons`` CLI flag sets this threshold to ``1/n``. When
      ``pool_polymorphisms`` is true, polymorphic sites in the outgroup are
      pooled into the ingroup counts.
    """

    dn: int
    """Non-synonymous fixed differences."""
    ds: int
    """Synonymous fixed differences."""
    pn: int
    """Non-synonymous polymorphisms in the ingroup (after frequency filtering)."""
    ps: int
    """Synonymous polymorphisms in the ingroup (after frequency filtering)."""
    p_value: float
    """Fisher's exact test p-value on the 2x2 contingency table."""
    ni: float | None
    """Neutrality Index ``(Pn/Ps) / (Dn/Ds)``; ``None`` if any count is zero."""
    alpha: float | None
    """Proportion of adaptive substitutions ``1 - NI``."""
    dos: float | None
    """Direction of Selection ``Dn/(Dn+Ds) - Pn/(Pn+Ps)`` (Stoletzki & Eyre-Walker 2011)."""
    ln: float | None = None
    """Total non-synonymous sites (Nei-Gojobori) over analyzed codons."""
    ls: float | None = None
    """Total synonymous sites (Nei-Gojobori) over analyzed codons."""
    omega: float | None = None
    """``(Dn/Ds) * (Ls/Ln)`` — dN/dS ratio. ``None`` when undefined.

    Note: this result class deliberately omits ``omega_a`` / ``omega_na``;
    the per-gene Smith & Eyre-Walker alpha is too noisy to give a useful
    rate decomposition. See ``docs/omega.rst`` for the rationale.
    """

    def __str__(self) -> str:
        ni_str = f"{self.ni:.4f}" if self.ni is not None else "NA"
        alpha_str = f"{self.alpha:.4f}" if self.alpha is not None else "NA"
        dos_str = f"{self.dos:.4f}" if self.dos is not None else "NA"
        omega_str = f"{self.omega:.4f}" if self.omega is not None else "NA"
        ls_str = f"{self.ls:.2f}" if self.ls is not None else "NA"
        ln_str = f"{self.ln:.2f}" if self.ln is not None else "NA"
        return (
            f"MK Test Results:\n"
            f"  Divergence:    Dn={self.dn}, Ds={self.ds}\n"
            f"  Polymorphism:  Pn={self.pn}, Ps={self.ps}\n"
            f"  Sites:         Ln={ln_str}, Ls={ls_str}\n"
            f"  Fisher's exact p-value: {self.p_value:.4g}\n"
            f"  Neutrality Index (NI):  {ni_str}\n"
            f"  Alpha (α):              {alpha_str}\n"
            f"  DoS:                    {dos_str}\n"
            f"  omega (dN/dS):          {omega_str}"
        )

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "dn": self.dn,
            "ds": self.ds,
            "pn": self.pn,
            "ps": self.ps,
            "p_value": self.p_value,
            "ni": self.ni,
            "alpha": self.alpha,
            "dos": self.dos,
            "ln": self.ln,
            "ls": self.ls,
            "omega": self.omega,
        }


def mk_test(
    ingroup: SequenceSet | str | Path,
    outgroup: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    pool_polymorphisms: bool = False,
    min_frequency: float = 0.0,
    no_singletons: bool = False,
) -> MKResult:
    """Perform the standard McDonald-Kreitman test.

    The MK test compares the ratio of non-synonymous to synonymous changes
    within species (polymorphism) to the ratio between species (divergence).

    Under neutrality, these ratios should be equal. Deviations suggest selection:
    - Excess divergence relative to polymorphism suggests positive selection
    - Excess polymorphism relative to divergence suggests segregating weakly deleterious polymorphisms

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup: SequenceSet or path to FASTA file for outgroup sequences
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        pool_polymorphisms: If True, count polymorphisms from both ingroup
            and outgroup (libsequence convention). If False (default), count
            only ingroup polymorphisms (DnaSP/original MK convention).
        min_frequency: Minimum derived allele frequency for a polymorphism
            to be counted (default 0.0, i.e., include all polymorphisms).
            Useful for filtering out rare variants that may be slightly
            deleterious.

    Returns:
        MKResult containing test statistics
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup, (str, Path)):
        outgroup = SequenceSet.from_fasta(outgroup, reading_frame, code)

    # Create aligned pair
    pair = AlignedPair(ingroup=ingroup, outgroup=outgroup, genetic_code=code)

    # Count divergence (fixed differences)
    dn = 0
    ds = 0
    processed_codons: set[int] = set()

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            if codon_idx in processed_codons:
                continue

            result = pair.classify_fixed_difference(codon_idx)
            if result is not None:
                nonsyn, syn = result
                dn += nonsyn
                ds += syn
                processed_codons.add(codon_idx)

    # Count polymorphisms
    pn = 0
    ps = 0

    if pool_polymorphisms:
        # Pooled mode: count polymorphisms from both populations (libsequence convention)
        poly_sites = pair.polymorphic_sites_pooled()
        classify_func = pair.classify_polymorphism_pooled
    else:
        # Standard mode: count only ingroup polymorphisms (DnaSP convention)
        poly_sites = pair.polymorphic_sites_ingroup()
        classify_func = pair.classify_polymorphism

    for codon_idx in poly_sites:
        # Skip if also a fixed difference (shouldn't happen, but be safe)
        if codon_idx in processed_codons:
            continue

        # Filtering needs the ancestral state, so a site whose state cannot be
        # determined is only dropped when a filter is actually in use.
        if min_frequency > 0 or no_singletons:
            state = ingroup.derived_state(outgroup, codon_idx)
            if state is None:
                continue
            derived_freq, derived_count = state
            if no_singletons and derived_count == 1:
                continue
            if derived_freq < min_frequency:
                continue

        result = classify_func(codon_idx)
        if result is not None:
            nonsyn, syn = result
            pn += nonsyn
            ps += syn

    # Calculate statistics
    p_val = fishers_exact(dn, ds, pn, ps)
    ni = neutrality_index(dn, ds, pn, ps)
    a = alpha(dn, ds, pn, ps)
    d = dos(dn, ds, pn, ps)
    ln, ls = pair.count_total_sites()
    omega, _, _ = omega_decomposition(dn, ds, ln, ls, a)

    return MKResult(
        dn=dn,
        ds=ds,
        pn=pn,
        ps=ps,
        p_value=p_val,
        ni=ni,
        alpha=a,
        dos=d,
        ln=ln,
        ls=ls,
        omega=omega,
    )


def mk_test_from_counts(
    dn: int,
    ds: int,
    pn: int,
    ps: int,
    ln: float | None = None,
    ls: float | None = None,
) -> MKResult:
    """Create MK test results from pre-computed counts.

    Useful for testing or when counts are already available.

    Args:
        dn: Non-synonymous divergence
        ds: Synonymous divergence
        pn: Non-synonymous polymorphisms
        ps: Synonymous polymorphisms
        ln: Total non-synonymous sites (optional, enables omega computation)
        ls: Total synonymous sites (optional, enables omega computation)

    Returns:
        MKResult with calculated statistics
    """
    p_val = fishers_exact(dn, ds, pn, ps)
    ni = neutrality_index(dn, ds, pn, ps)
    a = alpha(dn, ds, pn, ps)
    d = dos(dn, ds, pn, ps)
    omega, _, _ = omega_decomposition(dn, ds, ln, ls, a)

    return MKResult(
        dn=dn,
        ds=ds,
        pn=pn,
        ps=ps,
        p_value=p_val,
        ni=ni,
        alpha=a,
        dos=d,
        ln=ln,
        ls=ls,
        omega=omega,
    )
