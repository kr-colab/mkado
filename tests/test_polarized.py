"""Tests for the lineage counts of the polarized MK test."""

import pytest

from mkado.analysis.polarized import PolarizedMKResult, polarized_mk_test
from mkado.core.codons import GeneticCode
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


def run(ingroup: list[str], outgroup1: list[str], outgroup2: list[str], **kwargs):
    return polarized_mk_test(
        sequence_set(ingroup), sequence_set(outgroup1), sequence_set(outgroup2), **kwargs
    )


class TestUnpolarizedSites:
    """Sites the second outgroup cannot polarize are filtered and counted as in mk_test."""

    SINGLETON = ["AAA", "AAG", "AAA", "AAA"]

    def test_unpolarized_polymorphism_is_filtered(self) -> None:
        """The singleton AAG is removed by each filter, also when outgroup2 has no data."""
        assert polymorphism(run(self.SINGLETON, ["AAA"], ["NNN"])) == (0, 0, 0, 1)
        assert polymorphism(run(self.SINGLETON, ["AAA"], ["NNN"], no_singletons=True)) == (
            0,
            0,
            0,
            0,
        )
        assert polymorphism(run(self.SINGLETON, ["AAA"], ["NNN"], min_frequency=0.5)) == (
            0,
            0,
            0,
            0,
        )

    def test_stop_blocked_pair_counted_once(self) -> None:
        """Under table 2 every ordering from AAA to TGG passes through a stop codon."""
        result = run(["AAA"], ["TGG"], ["TGG"], genetic_code=GeneticCode(table_id=2))
        assert result.stop_blocked_pairs == 1

    def test_filtered_site_with_no_state_is_counted(self) -> None:
        """The ancestral AAC is in no ingroup copy, so the filter cannot judge the site."""
        result = run(["AAA", "AAG"], ["AAC"], ["AAC"], min_frequency=0.1)
        assert polymorphism(result) == (0, 0, 0, 0)
        assert result.unpolarized_sites == 1

    def test_each_mutation_filtered_on_its_own(self) -> None:
        """AAG is at 3/8 and passes; AGA is at 1/8 and does not."""
        result = run(["AAA"] * 4 + ["AAG"] * 3 + ["AGA"], ["AAA"], ["AAA"], min_frequency=0.2)
        assert polymorphism(result) == (0, 1, 0, 0)

    def test_outgroup2_polarizes_by_position(self) -> None:
        """Outgroup2 CAA shares the A at the polymorphic third position, though not the codon."""
        result = run(["AAA"] * 3 + ["AAG"], ["CAA"], ["CAA"], min_frequency=0.2)
        assert polymorphism(result) == (0, 1, 0, 0)

    def test_no_clean_outgroup1_codon_drops_the_polymorphism(self) -> None:
        """The codon is out of divergence and Ln/Ls, so its polymorphism is dropped and counted."""
        result = run(["AAA", "AAG"], ["NNN"], ["AAA"])
        assert polymorphism(result) == (0, 0, 0, 0)
        assert result.unpolarized_sites == 1

    def test_pooled_mode_uses_the_rule(self) -> None:
        result = run(["AAA", "AAG"], ["NNN"], ["AAA"], pool_polymorphisms=True)
        assert polymorphism(result) == (0, 0, 0, 0)
        assert result.unpolarized_sites == 1

    @pytest.mark.parametrize("outgroup2", [["AAG"], ["AAA", "AAG"]])
    def test_lineage_and_derived_state_use_one_ancestral_codon(self, outgroup2) -> None:
        """Outgroups 1 and 2 agree on AAG, so G is ancestral and the A at 6/8 is derived.

        A second outgroup2 allele must not change the ancestral base the filter reads.
        """
        result = run(["AAA"] * 6 + ["AAG", "AAC"], ["AAG"], outgroup2, min_frequency=0.2)
        assert polymorphism(result) == (0, 1, 0, 0)

    @pytest.mark.parametrize(
        ("ingroup", "outgroup1", "outgroup2", "expected"),
        [
            (SINGLETON, ["AAA"], ["NNN"], (0, 0, 0, 1)),
            (["AAA", "AAG"], ["AAC"], ["AAC"], (0, 1, 0, 0)),
            (["AAA"] * 4 + ["AAG"] * 3 + ["AGA"], ["AAA"], ["AAA"], (1, 1, 0, 0)),
            (["AAA"] * 3 + ["AAG"], ["CAA"], ["CAA"], (0, 1, 0, 0)),
            # Three alleles and an outgroup1 base shared with none: no state, but counted.
            (["AAA", "AAG", "AAC", "AAA"], ["AAT"], ["NNN"], (0, 0, 1, 1)),
        ],
    )
    def test_counts_without_a_filter_do_not_change(
        self, ingroup, outgroup1, outgroup2, expected
    ) -> None:
        result = run(ingroup, outgroup1, outgroup2)
        assert polymorphism(result) == expected
        assert result.unpolarized_sites == 0
