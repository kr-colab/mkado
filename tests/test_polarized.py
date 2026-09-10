"""Tests for the lineage counts of the polarized MK test."""

from mkado.analysis.polarized import PolarizedMKResult, polarized_mk_test
from tests.builders import sequence_set


def divergence(result: PolarizedMKResult) -> tuple[int, int, int, int, int, int]:
    """(Dn, Ds) on the ingroup lineage, then the outgroup lineage, then unpolarized."""
    return (
        result.dn_ingroup,
        result.ds_ingroup,
        result.dn_outgroup,
        result.ds_outgroup,
        result.dn_unpolarized,
        result.ds_unpolarized,
    )


def polymorphism(result: PolarizedMKResult) -> tuple[int, int, int, int]:
    """(Pn, Ps) on the ingroup lineage, then unpolarized."""
    return (result.pn_ingroup, result.ps_ingroup, result.pn_unpolarized, result.ps_unpolarized)


class TestFixedDifferences:
    """Each fixed difference lands on one lineage or in the unpolarized bucket."""

    def test_each_lineage_and_the_unpolarized_bucket(self) -> None:
        # Codon 0: AAC (Asn) in the ingroup, AAA (Lys) in both outgroups. A
        # replacement on the ingroup lineage.
        # Codon 1: AAA in the ingroup and outgroup2, AAG in outgroup1. A
        # synonymous change on the outgroup1 lineage.
        # Codon 2: AAA, AAG, and AAT across the three groups. Unpolarizable, so
        # it counts as a synonymous difference between the ingroup and outgroup1.
        # Codon 3: ATG everywhere.
        result = polarized_mk_test(
            ingroup=sequence_set(["AACAAAAAAATG", "AACAAAAAAATG"]),
            outgroup1=sequence_set(["AAAAAGAAGATG"]),
            outgroup2=sequence_set(["AAAAAAAATATG"]),
        )
        assert divergence(result) == (1, 0, 0, 1, 0, 1)
        assert polymorphism(result) == (0, 0, 0, 0)

    def test_polymorphic_second_outgroup_is_read_by_membership(self) -> None:
        # Codon 0: ingroup AAG, outgroup1 AAA, outgroup2 {AAA, AAT}. Outgroup2
        # carries the outgroup1 codon and not the ingroup's, so the ingroup changed.
        # Codon 1: ingroup AAG, outgroup1 AAA, outgroup2 {AAA, AAG}. Outgroup2
        # carries both, so the lineage cannot be chosen.
        result = polarized_mk_test(
            ingroup=sequence_set(["AAGAAG", "AAGAAG"]),
            outgroup1=sequence_set(["AAAAAA"]),
            outgroup2=sequence_set(["AAAAAA", "AATAAG"]),
        )
        assert divergence(result) == (0, 1, 0, 0, 0, 1)


class TestPolymorphisms:
    """Ingroup polymorphisms land on the ingroup lineage or in the unpolarized bucket."""

    INGROUP = ["AAA", "AAA", "AAG"]

    def test_polarized_polymorphism_counts_on_the_ingroup_lineage(self) -> None:
        result = polarized_mk_test(
            ingroup=sequence_set(self.INGROUP),
            outgroup1=sequence_set(["AAA"]),
            outgroup2=sequence_set(["AAA"]),
        )
        assert polymorphism(result) == (0, 1, 0, 0)

    def test_unpolarizable_polymorphism_goes_to_the_unpolarized_bucket(self) -> None:
        # Neither outgroup shares a codon with the ingroup or with the other.
        result = polarized_mk_test(
            ingroup=sequence_set(self.INGROUP),
            outgroup1=sequence_set(["AAC"]),
            outgroup2=sequence_set(["AAT"]),
        )
        assert polymorphism(result) == (0, 0, 0, 1)
        assert divergence(result) == (0, 0, 0, 0, 0, 0)

    def test_ancestral_polymorphism_goes_to_the_unpolarized_bucket(self) -> None:
        # Outgroup2 carries both ingroup alleles, so the polymorphism predates
        # the split and is not attributed to the ingroup lineage.
        result = polarized_mk_test(
            ingroup=sequence_set(self.INGROUP),
            outgroup1=sequence_set(["AAA"]),
            outgroup2=sequence_set(["AAA", "AAG"]),
        )
        assert polymorphism(result) == (0, 0, 0, 1)
