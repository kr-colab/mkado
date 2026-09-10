"""Codon classification over an aligned pair."""

from __future__ import annotations

from mkado.core.alignment import AlignedPair, PolarizedAlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import Sequence, SequenceSet

VERTEBRATE_MITO = GeneticCode(table_id=2)


def sequence_set(codons: list[str]) -> SequenceSet:
    return SequenceSet(sequences=[Sequence(name=f"s{i}", sequence=c) for i, c in enumerate(codons)])


def pair(ingroup: list[str], outgroup: list[str]) -> AlignedPair:
    return AlignedPair(
        ingroup=sequence_set(ingroup),
        outgroup=sequence_set(outgroup),
        genetic_code=DEFAULT_CODE,
    )


def pair_mito(ingroup: list[str], outgroup: list[str]) -> AlignedPair:
    return AlignedPair(
        ingroup=sequence_set(ingroup),
        outgroup=sequence_set(outgroup),
        genetic_code=VERTEBRATE_MITO,
    )


class TestClassifyPolymorphism:
    def test_two_codons(self):
        # AAA (Lys) to AAC (Asn) is one nonsynonymous change.
        assert pair(["AAA", "AAC"], ["AAA"]).classify_polymorphism(0) == (1, 0)

    def test_major_codon_is_the_most_common(self):
        """Every other codon is compared against the majority one."""
        ingroup = ["ACA"] * 5 + ["AAA", "AAC"]
        assert pair(ingroup, ["ACA"]).classify_polymorphism(0) == (1, 1)

    def test_monomorphic_site_is_not_a_polymorphism(self):
        assert pair(["AAA", "AAA"], ["AAA"]).classify_polymorphism(0) is None


class TestClassifyPolymorphismPooled:
    def test_two_codons_across_the_two_populations(self):
        assert pair(["AAA"], ["AAC"]).classify_polymorphism_pooled(0) == (1, 0)

    def test_major_codon_weighs_every_sequence_equally(self):
        """A large population must not be outvoted by a small one.

        Two ingroup sequences carry AAA while forty outgroup sequences split between
        ACA and AAC, so ACA is the majority codon of the pooled sample and each other
        codon is one step from it.
        """
        result = pair(["AAA"] * 2, ["ACA"] * 21 + ["AAC"] * 19)
        assert result.classify_polymorphism_pooled(0) == (1, 1)

    def test_agrees_with_classifying_the_concatenated_sample(self):
        """Pooling two populations is classifying the sequences of both together."""
        ingroup, outgroup = ["AAA"] * 2, ["ACA"] * 21 + ["AAC"] * 19
        pooled = pair(ingroup, outgroup).classify_polymorphism_pooled(0)
        together = pair(ingroup + outgroup, ["ACA"]).classify_polymorphism(0)
        assert pooled == together


class TestMultiAllelicCodons:
    """A codon carrying more than two alleles holds more than one mutation."""

    def test_two_alleles_at_one_position_are_two_mutations(self):
        """AAG is silent against AAA and AAC is a replacement, both at position 2."""
        ingroup = ["AAA"] * 3 + ["AAG", "AAC"]
        assert pair(ingroup, ["AAA"]).classify_polymorphism(0) == (1, 1)

    def test_answer_does_not_depend_on_sequence_order(self):
        """ACC and AAC both end in C at position 2, reached by different routes.

        Which route classifies that shared mutation depends on which allele is read
        first, so the order the alleles are visited has to come from the counts.
        """
        ingroup = ["AAA"] * 5 + ["ACC", "ACC", "AAC"]
        assert pair(ingroup, ["AAA"]).classify_polymorphism(0) == (1, 1)
        assert pair(list(reversed(ingroup)), ["AAA"]).classify_polymorphism(0) == (1, 1)

    def test_answer_does_not_depend_on_order_when_the_majority_ties(self):
        """Two codons tie for most common, so the tie cannot be broken by file order."""
        ingroup = ["AAA", "AAA", "AAC", "AAC", "AAG"]
        assert pair(ingroup, ["AAA"]).classify_polymorphism(0) == (1, 1)
        assert pair(list(reversed(ingroup)), ["AAA"]).classify_polymorphism(0) == (1, 1)

    def test_a_mutation_shared_by_two_alleles_counts_once(self):
        """ACC and ACG both carry the position 1 change, which happened once.

        Position 2 then carries two different bases, so the codon holds three
        mutations in total rather than four.
        """
        assert pair(["AAA"] * 3 + ["ACC", "ACG"], ["AAA"]).classify_polymorphism(0) == (1, 2)


class TestStopBlockedPairs:
    """Issue #90: a codon pair with no stop-free ordering is dropped, but counted.

    Under the vertebrate mitochondrial code (table 2), AAA/AAG (Lys) and
    TGA/TGG (Trp) differ at all three positions, and every one of the six
    orderings between them passes through AGA, AGG, or TAA -- all stops in
    that table. The pair is still real (both are ordinary sense codons), so
    it is worth counting even though it stays dropped from Dn/Ds/Pn/Ps.
    """

    def test_new_pair_starts_at_zero(self):
        assert pair_mito(["AAA"], ["TGG"]).stop_blocked_pairs == 0

    def test_classify_fixed_difference_counts_the_drop(self):
        result = pair_mito(["AAA"], ["TGG"])
        assert result.classify_fixed_difference(0) is None
        assert result.stop_blocked_pairs == 1

    def test_identical_codons_do_not_count(self):
        result = pair_mito(["AAA"], ["AAA"])
        assert result.classify_fixed_difference(0) is None
        assert result.stop_blocked_pairs == 0

    def test_missing_clean_codon_does_not_count(self):
        result = pair_mito(["NNN"], ["TGG"])
        assert result.classify_fixed_difference(0) is None
        assert result.stop_blocked_pairs == 0

    def test_counter_accumulates_across_codons(self):
        # One individual per group, two codons each: AAA/TGG and AAG/TGA are
        # both blocked pairs under table 2.
        result = AlignedPair(
            ingroup=SequenceSet(sequences=[Sequence(name="i1", sequence="AAAAAG")]),
            outgroup=SequenceSet(sequences=[Sequence(name="o1", sequence="TGGTGA")]),
            genetic_code=VERTEBRATE_MITO,
        )
        result.classify_fixed_difference(0)
        result.classify_fixed_difference(1)
        assert result.stop_blocked_pairs == 2

    def test_classify_against_major_single_other_counts_the_drop(self):
        result = pair_mito(["AAA", "TGG"], ["AAA"])
        assert result.classify_polymorphism(0) is None
        assert result.stop_blocked_pairs == 1

    def test_classify_against_major_multi_other_counts_only_the_blocked_one(self):
        # Major is AAA (Lys); TGG (Trp) is blocked, AAC (Asn) is a plain
        # single-position replacement.
        ingroup = ["AAA"] * 3 + ["TGG", "AAC"]
        result = pair_mito(ingroup, ["AAA"])
        nonsyn, syn = result.classify_polymorphism(0)
        assert (nonsyn, syn) == (1, 0)  # only AAC's replacement is counted
        assert result.stop_blocked_pairs == 1


class TestPolarizedStopBlockedPairs:
    """The same counter on PolarizedAlignedPair.polarize_fixed_difference."""

    def _pair(self, ingroup: str, outgroup1: str, outgroup2: str) -> PolarizedAlignedPair:
        return PolarizedAlignedPair(
            ingroup=sequence_set([ingroup]),
            outgroup=sequence_set([outgroup1]),
            outgroup2=sequence_set([outgroup2]),
            genetic_code=VERTEBRATE_MITO,
        )

    def test_blocked_polarized_difference_counts_the_drop(self):
        # Outgroup1 and outgroup2 agree on TGG, so it is ancestral and the
        # ingroup's AAA is derived -- a blocked pair on the ingroup lineage.
        pair = self._pair(ingroup="AAA", outgroup1="TGG", outgroup2="TGG")
        assert pair.polarize_fixed_difference(0) is None
        assert pair.stop_blocked_pairs == 1
