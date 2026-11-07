"""Alignment comparison utilities for MK test."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import (
    AlignmentLengthError,
    SequenceSet,
    ancestral_base,
    position_bases,
)

if TYPE_CHECKING:
    pass


def _count_path(path: list[tuple]) -> tuple[int, int]:
    """Split a mutational path, or a mutation list, into (replacement, synonymous) counts.

    Only the first item of each step, ``"R"`` or ``"S"``, is read.
    """
    nonsyn = sum(1 for step in path if step[0] == "R")
    return (nonsyn, len(path) - nonsyn)


def passes_frequency_filter(
    state: tuple[float, int] | None, min_frequency: float, no_singletons: bool
) -> bool:
    """Decide whether a polymorphism's (frequency, count) state clears the filters.

    ``AlignedPair.polarized_mutations``, ``SequenceSet.derived_state`` and
    ``AlignedPair.minor_allele_state`` all give this same (frequency, count)
    shape, so every counting loop can share this one accept/reject rule.
    """
    if state is None:
        return False
    freq, count = state
    if no_singletons and count == 1:
        return False
    return freq >= min_frequency


@dataclass
class AlignedPair:
    """Represents a pair of aligned sequence sets (ingroup and outgroup)."""

    ingroup: SequenceSet
    outgroup: SequenceSet
    genetic_code: GeneticCode = field(default_factory=lambda: DEFAULT_CODE)

    stop_blocked_pairs: int = field(default=0, init=False, repr=False)
    """Codon pairs where every mutational ordering passed through a stop codon
    and were dropped from the counts. See GeneticCode.get_path."""

    _internal_stops: set[int] = field(default_factory=set, init=False, repr=False)

    def __post_init__(self) -> None:
        self._check_aligned("outgroup", self.outgroup)
        # A stop codon is missing data: it adds nothing to divergence, polymorphism or the
        # site counts. The pair's code is used, because a set can carry the standard code
        # while the analysis uses another one.
        self.ingroup = self._masked(self.ingroup)
        self.outgroup = self._masked(self.outgroup)

    def _masked(self, group: SequenceSet) -> SequenceSet:
        """``group`` with its stop codons as NNN, recording the internal ones."""
        masked, internal = group.with_stops_masked(self.genetic_code)
        self._internal_stops |= internal
        return masked

    @property
    def trailing_bases(self) -> int:
        """Bases after the last full codon; every group has the ingroup's length and frame."""
        return self.ingroup.trailing_bases

    @property
    def internal_stop_codons(self) -> int:
        """Codons with an internal stop codon, one that is not its sequence's terminal stop."""
        return len(self._internal_stops)

    def _check_aligned(self, label: str, other: SequenceSet) -> None:
        """Raise if ``other`` and the ingroup differ in length or in reading frame.

        Bases are compared, not codons: a 6 bp and a 7 bp alignment both have two
        codons, but the extra base shows that they are not aligned. Each set reads its
        codons in its own frame, so sets in different frames would compare codons at
        different columns.
        """
        if self.ingroup.alignment_length != other.alignment_length:
            raise AlignmentLengthError(
                f"Alignment length mismatch: ingroup is {self.ingroup.alignment_length} bp, "
                f"{label} is {other.alignment_length} bp"
            )
        if self.ingroup.reading_frame != other.reading_frame:
            raise ValueError(
                f"Reading frame mismatch: ingroup is in frame {self.ingroup.reading_frame}, "
                f"{label} is in frame {other.reading_frame}"
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

    def minor_allele_state(self, codon_index: int) -> tuple[float, int] | None:
        """Frequency and copy count of the minor allele(s), pooled across both groups.

        Pooled mode has no ancestral state to compute a derived frequency from, so
        filtering there folds every non-major codon together instead. Ties for the
        major codon break alphabetically, the same as ``_classify_against_major``,
        so the result does not depend on sequence order in the file.

        Args:
            codon_index: Zero-based codon index.

        Returns:
            Tuple of (minor frequency, minor count), or None when fewer than two
            alleles are present.
        """
        counts = self.combined_codon_counts(codon_index)
        if len(counts) < 2:
            return None
        total = sum(counts.values())
        major_codon = min(counts, key=lambda c: (-counts[c], c))
        minor_count = total - counts[major_codon]
        return minor_count / total, minor_count

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

    def _get_path(self, codon1: str, codon2: str) -> list[tuple[str, int]] | None:
        """Look up a mutational path, counting and reporting a stop-blocked pair as None."""
        path = self.genetic_code.get_path(codon1, codon2)
        if not path:
            self.stop_blocked_pairs += 1
            return None
        return path

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

        path = self._get_path(in_codon, out_codon)
        if path is None:
            return None

        return _count_path(path)

    def _mutations_against_major(self, counts: dict[str, int]) -> list[tuple[str, int, str, str]]:
        """List a codon's mutations against the most common allele.

        Each mutation is ``(change_type, position, base, major_base)``: the base this
        mutation gives at ``position``, and the major allele's base there.

        Every allele is compared with the majority codon. A mutation is identified by
        the position it changes and the base it produces, so a change two alleles share
        is counted once while two different bases at one position count as two.

        Alleles are ordered most common first, ties broken alphabetically, and the first
        is the one the rest are compared against. Two alleles can reach the same base at
        one position through different intermediate codons, which classify differently,
        so a fixed order is what keeps the result independent of the order of sequences
        in the file.

        ``change_type`` is ``"R"`` for a replacement and ``"S"`` for a synonymous
        change, and ``position`` is the codon position (0, 1 or 2) it changes.

        Args:
            counts: Codon counts at this position over the sample being classified.

        Returns:
            The mutations, empty when fewer than two alleles are present.
        """
        if len(counts) < 2:
            return []

        major_codon, *others = sorted(counts, key=lambda c: (-counts[c], c))

        mutations: list[tuple[str, int, str, str]] = []
        counted: set[tuple[int, str]] = set()
        for codon in others:
            path = self._get_path(major_codon, codon)
            if path is None:
                # No stop-free ordering exists; this allele is dropped, not the site.
                continue
            for change_type, pos in path:
                mutation = (pos, codon[pos])
                if mutation in counted:
                    continue
                counted.add(mutation)
                mutations.append((change_type, pos, codon[pos], major_codon[pos]))
        return mutations

    def _classify_against_major(self, counts: dict[str, int]) -> tuple[int, int] | None:
        """Count a codon's mutations as (non_synonymous, synonymous), or None if not polymorphic.

        See ``_mutations_against_major`` for how the mutations are found.
        """
        if len(counts) < 2:
            return None
        return _count_path(self._mutations_against_major(counts))

    def polarized_mutations(
        self,
        codon_index: int,
        pooled: bool = False,
        ancestral_codons: set[str] | None = None,
    ) -> tuple[list[tuple[str, tuple[float, int] | None]], int]:
        """Each polymorphic mutation at a codon with its derived state, and how many were dropped.

        A mutation is ``(change_type, state)``: ``change_type`` is ``"R"`` or ``"S"``, and
        ``state`` is the (frequency, count) of the mutation's derived base at its position,
        or None when the outgroup shares no base there. The ancestral base is chosen as in
        ``SequenceSet.position_derived_state``, but at a position with three bases each
        derived base has its own frequency. Polymorphism is counted only where both groups have a clean codon, the same
        rule as for divergence and ``count_total_sites``; elsewhere every mutation is
        dropped and returned as the drop count. Pooled mode lists the mutations over both
        groups, as ``classify_polymorphism_pooled`` does, but polarizes the ingroup.

        ``ancestral_codons`` gives the codons whose bases name the ancestral base, the
        outgroup's clean codons by default. The drop rule always reads the pair's own
        outgroup, because that is the rule for divergence and the site counts.
        """
        ingroup_counts, total = self.ingroup.site_codon_counts(codon_index)
        counts = self.combined_codon_counts(codon_index) if pooled else ingroup_counts
        mutations = self._mutations_against_major(counts)
        if not (ingroup_counts and self.outgroup.codon_set_clean(codon_index)):
            return [], len(mutations)
        source_codons = (
            self.outgroup.codon_set_clean(codon_index)
            if ancestral_codons is None
            else ancestral_codons
        )
        by_position = {}
        for position in {position for _, position, _, _ in mutations}:
            bases = position_bases(ingroup_counts, position)
            by_position[position] = (bases, ancestral_base(bases, source_codons, position))
        polarized = []
        for change_type, position, base, major_base in mutations:
            bases, ancestral = by_position[position]
            state = None
            if ancestral is not None:
                # When the mutation's base is the ancestral one, the major allele carries
                # the derived base.
                derived = major_base if base == ancestral else base
                count = bases.get(derived, 0)
                state = (count / total, count)
            polarized.append((change_type, state))
        return polarized, 0

    def filtered_polymorphism(
        self,
        codon_index: int,
        min_frequency: float = 0.0,
        no_singletons: bool = False,
        *,
        pooled: bool = False,
        ancestral_codons: set[str] | None = None,
    ) -> tuple[int, int, int]:
        """(nonsyn, syn, dropped) for the polymorphism at a codon, after the frequency filters.

        Each mutation is filtered by its own position's derived state. Without a filter no
        ancestral state is needed, so a third-allele site still counts; with one, a mutation
        with no derived state cannot be judged and is dropped. Pooled mode has no ancestral
        state, so it folds both groups together and filters the codon on its minor allele.
        ``dropped`` also counts the mutations of a codon with no clean outgroup codon.
        """
        mutations, dropped = self.polarized_mutations(
            codon_index, pooled=pooled, ancestral_codons=ancestral_codons
        )
        if min_frequency > 0 or no_singletons:
            if pooled:
                state = self.minor_allele_state(codon_index)
                if not passes_frequency_filter(state, min_frequency, no_singletons):
                    return 0, 0, dropped
            else:
                dropped += sum(1 for _, state in mutations if state is None)
                mutations = [
                    m
                    for m in mutations
                    if passes_frequency_filter(m[1], min_frequency, no_singletons)
                ]
        nonsyn, syn = _count_path(mutations)
        return nonsyn, syn, dropped

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

        Follows the libsequence convention of drawing codons from both the ingroup
        and the outgroup. Their counts are pooled, so a codon's multiplicity decides
        which one the others are compared against.

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

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.outgroup2 is not None:
            self._check_aligned("outgroup2", self.outgroup2)
            self.outgroup2 = self._masked(self.outgroup2)

    def polarize_fixed_difference(self, codon_index: int) -> tuple[str, tuple[int, int]] | None:
        """Polarize a fixed difference to determine which lineage changed.

        The second outgroup names the ancestral codon. It may be polymorphic at
        the site, so it is read by membership: if it carries the first
        outgroup's codon, the ingroup changed; if it carries the ingroup's
        codon, the first outgroup changed; if it carries both or neither, the
        lineage cannot be chosen.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (lineage, (nonsyn, syn)) where lineage is 'ingroup' or
            'outgroup', or None if cannot be polarized. A stop-blocked pair also gives
            None and is counted here, so a caller must not classify it again.
        """
        if self.outgroup2 is None or not self.is_fixed_between(codon_index):
            return None

        (in_codon,) = self.ingroup.codon_set_clean(codon_index)
        (out_codon,) = self.outgroup.codon_set_clean(codon_index)
        out2_codons = self.outgroup2.codon_set_clean(codon_index)

        out_is_ancestral = out_codon in out2_codons
        in_is_ancestral = in_codon in out2_codons
        if out_is_ancestral == in_is_ancestral:
            return None

        if out_is_ancestral:
            ancestral, derived, lineage = out_codon, in_codon, "ingroup"
        else:
            ancestral, derived, lineage = in_codon, out_codon, "outgroup"

        path = self._get_path(ancestral, derived)
        if path is None:
            return None

        return (lineage, _count_path(path))

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
        if self.ingroup_ancestral_codons(codon_index) is None:
            return None
        return self.classify_polymorphism(codon_index)

    def ingroup_ancestral_codons(self, codon_index: int) -> set[str] | None:
        """The ancestral codons of an ingroup polymorphism on the ingroup lineage, or None.

        None when the polymorphism cannot be put on the ingroup lineage; see
        ``polarize_ingroup_polymorphism`` for the rules. The returned codons are the ones
        that name the ancestral base at each position when the polymorphism is filtered.
        """
        if self.outgroup2 is None:
            return None

        in_codons = self.ingroup.codon_set_clean(codon_index)
        out1_codons = self.outgroup.codon_set_clean(codon_index)
        out2_codons = self.outgroup2.codon_set_clean(codon_index)

        if len(in_codons) < 2 or not out2_codons:
            return None  # Not polymorphic in the ingroup, or no outgroup2 data
        if in_codons <= out2_codons:
            return None  # Ancestral polymorphism: not on the ingroup lineage

        # Preferred: outgroup1 and outgroup2 agree on the ancestral codon.
        if shared := out1_codons & out2_codons:
            return shared
        # Alternative: outgroup2 shares some (but not all) alleles with the ingroup.
        return (in_codons & out2_codons) or None
