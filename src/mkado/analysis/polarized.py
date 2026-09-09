"""Polarized McDonald-Kreitman test implementation."""

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
from mkado.core.alignment import PolarizedAlignedPair, passes_frequency_filter
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class PolarizedMKResult:
    """Results from a polarized McDonald-Kreitman test.

    The polarized test uses a second outgroup (``outgroup2``) to assign each
    fixed difference and each ingroup polymorphism to one of two lineages —
    the ingroup or the first outgroup — by majority rule on the inferred
    ancestral state. Sites where the ancestral state cannot be determined
    are reported under ``*_unpolarized``.

    Pn and Ps count derived alleles segregating in the ingroup (subject to
    ``min_frequency`` filtering, identical to the standard MK test) and
    are split by lineage assignment. In pooled mode there is no ancestral
    state to assign a lineage or a derived allele from, so filtering there
    uses the folded minor-allele frequency over the combined ingroup+outgroup
    sample instead, and lineage assignment does not apply.
    """

    dn_ingroup: int
    """Non-synonymous fixed differences assigned to the ingroup lineage."""
    ds_ingroup: int
    """Synonymous fixed differences assigned to the ingroup lineage."""
    pn_ingroup: int
    """Non-synonymous polymorphisms whose derived allele arose on the ingroup lineage."""
    ps_ingroup: int
    """Synonymous polymorphisms whose derived allele arose on the ingroup lineage."""

    dn_outgroup: int
    """Non-synonymous fixed differences assigned to the first outgroup lineage."""
    ds_outgroup: int
    """Synonymous fixed differences assigned to the first outgroup lineage."""

    dn_unpolarized: int
    """Non-synonymous fixed differences where lineage could not be determined."""
    ds_unpolarized: int
    """Synonymous fixed differences where lineage could not be determined."""
    pn_unpolarized: int
    """Non-synonymous polymorphisms where ancestral state could not be determined."""
    ps_unpolarized: int
    """Synonymous polymorphisms where ancestral state could not be determined."""

    p_value_ingroup: float
    """Fisher's exact test p-value on the ingroup-lineage 2x2 table."""
    ni_ingroup: float | None
    """Neutrality Index for the ingroup lineage."""
    alpha_ingroup: float | None
    """Proportion of adaptive substitutions on the ingroup lineage."""
    dos_ingroup: float | None
    """Direction of Selection for the ingroup lineage."""
    ln: float | None = None
    """Total non-synonymous sites (Nei-Gojobori) over analyzed codons."""
    ls: float | None = None
    """Total synonymous sites (Nei-Gojobori) over analyzed codons."""
    omega: float | None = None
    """``(Dn_ingroup/Ds_ingroup) * (Ls/Ln)`` — ingroup-lineage dN/dS.

    Note: this result class deliberately omits ``omega_a`` / ``omega_na``;
    per-gene polarized alpha is too noisy to give a useful rate
    decomposition. See ``docs/omega.rst`` for the rationale.
    """

    def __str__(self) -> str:
        ni_str = f"{self.ni_ingroup:.4f}" if self.ni_ingroup is not None else "NA"
        alpha_str = f"{self.alpha_ingroup:.4f}" if self.alpha_ingroup is not None else "NA"
        dos_str = f"{self.dos_ingroup:.4f}" if self.dos_ingroup is not None else "NA"
        omega_str = f"{self.omega:.4f}" if self.omega is not None else "NA"
        ls_str = f"{self.ls:.2f}" if self.ls is not None else "NA"
        ln_str = f"{self.ln:.2f}" if self.ln is not None else "NA"
        return (
            f"Polarized MK Test Results:\n"
            f"  Ingroup Lineage:\n"
            f"    Divergence:    Dn={self.dn_ingroup}, Ds={self.ds_ingroup}\n"
            f"    Polymorphism:  Pn={self.pn_ingroup}, Ps={self.ps_ingroup}\n"
            f"    Fisher's p-value: {self.p_value_ingroup:.4g}\n"
            f"    NI: {ni_str}, Alpha: {alpha_str}, DoS: {dos_str}\n"
            f"    omega (dN/dS): {omega_str}\n"
            f"  Sites: Ln={ln_str}, Ls={ls_str}\n"
            f"  Outgroup Lineage:\n"
            f"    Divergence:    Dn={self.dn_outgroup}, Ds={self.ds_outgroup}\n"
            f"  Unpolarized:\n"
            f"    Divergence:    Dn={self.dn_unpolarized}, Ds={self.ds_unpolarized}\n"
            f"    Polymorphism:  Pn={self.pn_unpolarized}, Ps={self.ps_unpolarized}"
        )

    def to_dict(self) -> dict:
        """Convert results to a dictionary."""
        return {
            "ingroup": {
                "dn": self.dn_ingroup,
                "ds": self.ds_ingroup,
                "pn": self.pn_ingroup,
                "ps": self.ps_ingroup,
                "p_value": self.p_value_ingroup,
                "ni": self.ni_ingroup,
                "alpha": self.alpha_ingroup,
                "dos": self.dos_ingroup,
                "omega": self.omega,
            },
            "outgroup": {
                "dn": self.dn_outgroup,
                "ds": self.ds_outgroup,
            },
            "unpolarized": {
                "dn": self.dn_unpolarized,
                "ds": self.ds_unpolarized,
                "pn": self.pn_unpolarized,
                "ps": self.ps_unpolarized,
            },
            "ln": self.ln,
            "ls": self.ls,
        }


def polarized_mk_test(
    ingroup: SequenceSet | str | Path,
    outgroup1: SequenceSet | str | Path,
    outgroup2: SequenceSet | str | Path,
    reading_frame: int = 1,
    genetic_code: GeneticCode | None = None,
    pool_polymorphisms: bool = False,
    min_frequency: float = 0.0,
    no_singletons: bool = False,
) -> PolarizedMKResult:
    """Perform a polarized McDonald-Kreitman test.

    Uses a second outgroup to determine the ancestral state and assign
    mutations to specific lineages.

    Args:
        ingroup: SequenceSet or path to FASTA file for ingroup sequences
        outgroup1: SequenceSet or path to FASTA file for first outgroup
        outgroup2: SequenceSet or path to FASTA file for second outgroup
            (used for polarization)
        reading_frame: Reading frame (1, 2, or 3)
        genetic_code: Genetic code for translation (uses standard if None)
        pool_polymorphisms: If True, count polymorphisms from both ingroup
            and outgroup1 (libsequence convention). If False (default), count
            only ingroup polymorphisms (DnaSP/original MK convention). Pooled
            mode has no ancestral state, so ``min_frequency``/``no_singletons``
            filter on the folded minor-allele frequency over the combined
            sample rather than a derived-allele frequency.
        min_frequency: Minimum derived (or, in pooled mode, minor) allele
            frequency for a polymorphism to be counted (default 0.0, i.e.,
            include all polymorphisms). Useful for filtering out rare variants
            that may be slightly deleterious.

    Returns:
        PolarizedMKResult containing test statistics
    """
    code = genetic_code or DEFAULT_CODE

    # Load sequences if paths provided
    if isinstance(ingroup, (str, Path)):
        ingroup = SequenceSet.from_fasta(ingroup, reading_frame, code)
    if isinstance(outgroup1, (str, Path)):
        outgroup1 = SequenceSet.from_fasta(outgroup1, reading_frame, code)
    if isinstance(outgroup2, (str, Path)):
        outgroup2 = SequenceSet.from_fasta(outgroup2, reading_frame, code)

    # Create polarized aligned pair
    pair = PolarizedAlignedPair(
        ingroup=ingroup,
        outgroup=outgroup1,
        outgroup2=outgroup2,
        genetic_code=code,
    )

    # Count divergence with polarization
    dn_in = 0
    ds_in = 0
    dn_out = 0
    ds_out = 0
    dn_unpol = 0
    ds_unpol = 0

    for codon_idx in range(pair.num_codons):
        if pair.is_fixed_between(codon_idx):
            polarization = pair.polarize_fixed_difference(codon_idx)

            if polarization is None:
                # Cannot polarize - count as unpolarized
                result = pair.classify_fixed_difference(codon_idx)
                if result is not None:
                    nonsyn, syn = result
                    dn_unpol += nonsyn
                    ds_unpol += syn
            else:
                lineage, (nonsyn, syn) = polarization
                if lineage == "ingroup":
                    dn_in += nonsyn
                    ds_in += syn
                else:
                    dn_out += nonsyn
                    ds_out += syn

    # Count polymorphisms with polarization
    pn_in = 0
    ps_in = 0
    pn_unpol = 0
    ps_unpol = 0

    if pool_polymorphisms:
        # Pooled mode: count polymorphisms from both populations (libsequence convention)
        # Note: polarization is not applied in pooled mode as it's unclear which
        # lineage to attribute shared polymorphisms to
        poly_sites = pair.polymorphic_sites_pooled()
        for codon_idx in poly_sites:
            # Pooled mode has no ancestral state, so filtering folds ingroup and
            # outgroup together and uses the minor allele instead of the derived one.
            if min_frequency > 0 or no_singletons:
                state = pair.minor_allele_state(codon_idx)
                if not passes_frequency_filter(state, min_frequency, no_singletons):
                    continue

            result = pair.classify_polymorphism_pooled(codon_idx)
            if result is not None:
                nonsyn, syn = result
                pn_in += nonsyn
                ps_in += syn
    else:
        # Standard mode: count only ingroup polymorphisms (DnaSP convention)
        # Apply polarization to determine if polymorphism arose on ingroup lineage
        poly_sites = pair.polymorphic_sites_ingroup()

        for codon_idx in poly_sites:
            # First, try to polarize the polymorphism
            result = pair.polarize_ingroup_polymorphism(codon_idx)

            if result is None:
                # Could not polarize - count as unpolarized
                unpol_result = pair.classify_polymorphism(codon_idx)
                if unpol_result is not None:
                    nonsyn, syn = unpol_result
                    pn_unpol += nonsyn
                    ps_unpol += syn
                continue

            # Polymorphism is polarized to ingroup lineage. The second outgroup carries
            # the ancestral state here, as it does for the polarization itself.
            if min_frequency > 0 or no_singletons:
                state = ingroup.derived_state(outgroup2, codon_idx)
                if not passes_frequency_filter(state, min_frequency, no_singletons):
                    continue

            nonsyn, syn = result
            pn_in += nonsyn
            ps_in += syn

    # Calculate statistics for ingroup lineage
    p_val = fishers_exact(dn_in, ds_in, pn_in, ps_in)
    ni = neutrality_index(dn_in, ds_in, pn_in, ps_in)
    a = alpha(dn_in, ds_in, pn_in, ps_in)
    d = dos(dn_in, ds_in, pn_in, ps_in)
    ln, ls = pair.count_total_sites()
    omega, _, _ = omega_decomposition(dn_in, ds_in, ln, ls, a)

    return PolarizedMKResult(
        dn_ingroup=dn_in,
        ds_ingroup=ds_in,
        pn_ingroup=pn_in,
        ps_ingroup=ps_in,
        dn_outgroup=dn_out,
        ds_outgroup=ds_out,
        dn_unpolarized=dn_unpol,
        ds_unpolarized=ds_unpol,
        pn_unpolarized=pn_unpol,
        ps_unpolarized=ps_unpol,
        p_value_ingroup=p_val,
        ni_ingroup=ni,
        alpha_ingroup=a,
        dos_ingroup=d,
        ln=ln,
        ls=ls,
        omega=omega,
    )
