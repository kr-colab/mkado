"""Codon classification over an aligned pair."""

from __future__ import annotations

from mkado.core.alignment import AlignedPair
from mkado.core.codons import DEFAULT_CODE
from mkado.core.sequences import Sequence, SequenceSet


def sequence_set(codons: list[str]) -> SequenceSet:
    return SequenceSet(sequences=[Sequence(name=f"s{i}", sequence=c) for i, c in enumerate(codons)])


def pair(ingroup: list[str], outgroup: list[str]) -> AlignedPair:
    return AlignedPair(
        ingroup=sequence_set(ingroup),
        outgroup=sequence_set(outgroup),
        genetic_code=DEFAULT_CODE,
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
