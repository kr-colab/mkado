"""Alignment comparison utilities for MK test."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import SequenceSet

if TYPE_CHECKING:
    pass


@dataclass
class AlignedPair:
    """Represents a pair of aligned sequence sets (ingroup and outgroup)."""

    ingroup: SequenceSet
    outgroup: SequenceSet
    genetic_code: GeneticCode = field(default_factory=lambda: DEFAULT_CODE)

    def __post_init__(self) -> None:
        if self.ingroup.num_codons != self.outgroup.num_codons:
            raise ValueError(
                f"Alignment length mismatch: ingroup has {self.ingroup.num_codons} codons, "
                f"outgroup has {self.outgroup.num_codons} codons"
            )

    @property
    def num_codons(self) -> int:
        """Number of codons in the alignment."""
        return self.ingroup.num_codons

    def combined_codon_set(self, codon_index: int) -> set[str]:
        """Get all unique codons at a position from both groups.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique codon strings from both ingroup and outgroup
        """
        return self.ingroup.codon_set(codon_index) | self.outgroup.codon_set(codon_index)

    def combined_codon_counts(self, codon_index: int) -> dict[str, int]:
        """Count each clean codon at a position over both groups together.

        Adding the two frequency spectra instead would weight the groups equally
        regardless of how many sequences each holds.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Dict mapping codon strings to the number of sequences carrying them
        """
        ingroup_counts, _ = self.ingroup.site_codon_counts(codon_index)
        outgroup_counts, _ = self.outgroup.site_codon_counts(codon_index)
        combined = dict(ingroup_counts)
        for codon, count in outgroup_counts.items():
            combined[codon] = combined.get(codon, 0) + count
        return combined

    def combined_codon_set_clean(self, codon_index: int) -> set[str]:
        """Get all clean unique codons at a position from both groups.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique valid codon strings
        """
        return self.ingroup.codon_set_clean(codon_index) | self.outgroup.codon_set_clean(
            codon_index
        )

    def is_fixed_between(self, codon_index: int) -> bool:
        """Check if a codon position is fixed differently between groups.

        A position is a fixed difference if:
        - The ingroup has a single codon (or all share the same)
        - The outgroup has a single codon (or all share the same)
        - The ingroup and outgroup codons are different

        Args:
            codon_index: Zero-based codon index

        Returns:
            True if this is a fixed difference between groups
        """
        in_codons = self.ingroup.codon_set_clean(codon_index)
        out_codons = self.outgroup.codon_set_clean(codon_index)

        if not in_codons or not out_codons:
            return False

        # Fixed if each group has one codon and they're different
        if len(in_codons) == 1 and len(out_codons) == 1:
            return in_codons != out_codons

        return False

    def is_polymorphic_within_ingroup(self, codon_index: int) -> bool:
        """Check if a codon position is polymorphic within the ingroup only.

        Args:
            codon_index: Zero-based codon index

        Returns:
            True if polymorphic within ingroup but not a fixed difference
        """
        return self.ingroup.is_polymorphic(codon_index)

    def fixed_differences(self) -> list[int]:
        """Get all codon indices with fixed differences between groups.

        Returns:
            List of codon indices
        """
        return [i for i in range(self.num_codons) if self.is_fixed_between(i)]

    def polymorphic_sites_ingroup(self) -> list[int]:
        """Get all codon indices polymorphic within the ingroup.

        Returns:
            List of codon indices
        """
        return self.ingroup.polymorphic_codons()

    def polymorphic_sites_outgroup(self) -> list[int]:
        """Get all codon indices polymorphic within the outgroup.

        Returns:
            List of codon indices
        """
        return self.outgroup.polymorphic_codons()

    def polymorphic_sites_pooled(self) -> list[int]:
        """Get all codon indices polymorphic in either population (union).

        This follows the libsequence convention of pooling polymorphisms
        from both populations.

        Returns:
            List of unique codon indices (sorted)
        """
        ingroup_poly = set(self.ingroup.polymorphic_codons())
        outgroup_poly = set(self.outgroup.polymorphic_codons())
        return sorted(ingroup_poly | outgroup_poly)

    def count_total_sites(self) -> tuple[float, float]:
        """Aggregate total non-synonymous (Ln) and synonymous (Ls) sites.

        For each codon position with at least one clean codon in each group,
        averages the Nei-Gojobori per-codon synonymous site counts within
        each group, then averages the two group means. Sums across positions
        and returns ``(Ln, Ls)`` with ``Ln = 3 * n_codons - Ls``.
        """
        n_codons = self.num_codons
        ingroup = self.ingroup
        outgroup = self.outgroup
        cs = self.genetic_code.count_synonymous_sites

        ls = 0.0
        used = 0
        for i in range(n_codons):
            in_codons_i = ingroup.codon_set_clean(i)
            out_codons_i = outgroup.codon_set_clean(i)
            if not in_codons_i or not out_codons_i:
                continue
            in_ls = sum(cs(c) for c in in_codons_i) / len(in_codons_i)
            out_ls = sum(cs(c) for c in out_codons_i) / len(out_codons_i)
            ls += (in_ls + out_ls) / 2.0
            used += 1

        ln = (3.0 * used) - ls
        return (ln, ls)

    def classify_fixed_difference(self, codon_index: int) -> tuple[int, int] | None:
        """Classify a fixed difference as synonymous/non-synonymous.

        Counts changes along the path with the fewest replacements, which is
        MKado's rule for a pair differing at more than one position. See
        :ref:`counting-differences`.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if
            not a valid fixed difference
        """
        in_codons = self.ingroup.codon_set_clean(codon_index)
        out_codons = self.outgroup.codon_set_clean(codon_index)

        if not in_codons or not out_codons:
            return None

        # Get representative codons
        in_codon = next(iter(in_codons))
        out_codon = next(iter(out_codons))

        if in_codon == out_codon:
            return None

        path = self.genetic_code.get_path(in_codon, out_codon)
        if not path:
            return None

        nonsyn = sum(1 for change_type, _ in path if change_type == "R")
        syn = sum(1 for change_type, _ in path if change_type == "S")

        return (nonsyn, syn)

    def _classify_against_major(self, counts: dict[str, int]) -> tuple[int, int] | None:
        """Classify a codon's alleles against the most common one.

        Every allele is compared with the majority codon, and each nucleotide position
        contributes at most once, so one mutation is not counted down two paths.

        Args:
            counts: Codon counts at this position over the sample being classified.

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if not a valid
            polymorphism.
        """
        codons = list(counts)
        if len(codons) < 2:
            return None

        if len(codons) == 2:
            path = self.genetic_code.get_path(codons[0], codons[1])
            if not path:
                return None
            nonsyn = sum(1 for change_type, _ in path if change_type == "R")
            syn = sum(1 for change_type, _ in path if change_type == "S")
            return (nonsyn, syn)

        total_nonsyn = 0
        total_syn = 0
        counted_positions: set[int] = set()
        major_codon = max(counts, key=lambda c: counts[c])
        for codon in codons:
            if codon == major_codon:
                continue
            for change_type, pos in self.genetic_code.get_path(major_codon, codon):
                if pos not in counted_positions:
                    if change_type == "R":
                        total_nonsyn += 1
                    else:
                        total_syn += 1
                    counted_positions.add(pos)

        return (total_nonsyn, total_syn)

    def classify_polymorphism(self, codon_index: int) -> tuple[int, int] | None:
        """Classify a polymorphism as synonymous/non-synonymous.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if
            not a valid polymorphism
        """
        counts, _ = self.ingroup.site_codon_counts(codon_index)
        return self._classify_against_major(counts)

    def classify_polymorphism_pooled(self, codon_index: int) -> tuple[int, int] | None:
        """Classify a polymorphism using codons from both populations.

        This follows the libsequence convention of using all unique codons
        from both ingroup and outgroup when classifying polymorphisms.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count), or None if
            not a valid polymorphism
        """
        return self._classify_against_major(self.combined_codon_counts(codon_index))


@dataclass
class PolarizedAlignedPair(AlignedPair):
    """Aligned pair with a second outgroup for polarization."""

    outgroup2: SequenceSet | None = None

    def polarize_fixed_difference(self, codon_index: int) -> tuple[str, tuple[int, int]] | None:
        """Polarize a fixed difference to determine which lineage changed.

        Uses the second outgroup to determine the ancestral state.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (lineage, (nonsyn, syn)) where lineage is 'ingroup' or
            'outgroup', or None if cannot be polarized
        """
        if self.outgroup2 is None:
            return None

        in_codons = self.ingroup.codon_set_clean(codon_index)
        out_codons = self.outgroup.codon_set_clean(codon_index)
        out2_codons = self.outgroup2.codon_set_clean(codon_index)

        if not in_codons or not out_codons or not out2_codons:
            return None

        in_codon = next(iter(in_codons))
        out_codon = next(iter(out_codons))
        out2_codon = next(iter(out2_codons))

        # Determine ancestral state
        if out_codon == out2_codon:
            # Outgroup agrees - ingroup changed
            ancestral = out_codon
            derived = in_codon
            lineage = "ingroup"
        elif in_codon == out2_codon:
            # Ingroup matches outgroup2 - outgroup1 changed
            ancestral = in_codon
            derived = out_codon
            lineage = "outgroup"
        else:
            # Cannot polarize
            return None

        if ancestral == derived:
            return None

        path = self.genetic_code.get_path(ancestral, derived)
        if not path:
            return None

        nonsyn = sum(1 for change_type, _ in path if change_type == "R")
        syn = sum(1 for change_type, _ in path if change_type == "S")

        return (lineage, (nonsyn, syn))

    def polarize_ingroup_polymorphism(self, codon_index: int) -> tuple[int, int] | None:
        """Polarize and classify an ingroup polymorphism.

        Uses outgroup2 to determine if the polymorphism arose on the ingroup
        lineage. Following the convention from mkTest.rb:

        1. If outgroup2 has ALL ingroup alleles → ancestral polymorphism,
           cannot attribute to ingroup lineage
        2. If outgroup2 shares allele with outgroup1 → the shared allele is
           ancestral, ingroup has derived allele(s) → ingroup polymorphism
        3. If outgroup2 shares some (but not all) alleles with ingroup →
           shared allele is ancestral → ingroup polymorphism
        4. Otherwise → cannot polarize

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (non_synonymous_count, synonymous_count) if the
            polymorphism is derived in ingroup, None if cannot polarize
        """
        if self.outgroup2 is None:
            return None

        in_codons = self.ingroup.codon_set_clean(codon_index)
        out1_codons = self.outgroup.codon_set_clean(codon_index)
        out2_codons = self.outgroup2.codon_set_clean(codon_index)

        if not in_codons or len(in_codons) < 2:
            return None  # Not polymorphic in ingroup

        if not out2_codons:
            return None  # No outgroup2 data - cannot polarize

        # Check if outgroup2 has all ingroup alleles (ancestral polymorphism)
        if in_codons <= out2_codons:
            return None  # Ancestral polymorphism - don't count for ingroup lineage

        # Check if we can determine ancestral state
        # Preferred: outgroup1 and outgroup2 agree (both have same allele)
        if out1_codons and (out1_codons & out2_codons):
            # outgroup1 and outgroup2 share an allele - that's the ancestral state
            return self.classify_polymorphism(codon_index)

        # Alternative: outgroup2 shares some (but not all) alleles with ingroup
        if in_codons & out2_codons:
            # The shared allele is ancestral, others are derived in ingroup
            return self.classify_polymorphism(codon_index)

        return None  # Cannot polarize
