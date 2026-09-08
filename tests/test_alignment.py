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
