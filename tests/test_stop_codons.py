"""Stop codons are missing data: they add nothing to divergence, polymorphism or site counts."""

from __future__ import annotations

import pytest

from mkado.analysis.asymptotic import extract_polymorphism_data
from mkado.analysis.mk_test import mk_test
from mkado.analysis.polarized import polarized_mk_test
from mkado.core.alignment import PolarizedAlignedPair
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from tests.builders import sequence_set

VERTEBRATE_MITO = GeneticCode(table_id=2)


def counts(result) -> tuple[int, int, int, int]:
    return result.dn, result.ds, result.pn, result.ps


@pytest.mark.parametrize(
    ("ingroup", "outgroup", "expected"),
    [
        # A difference between two stop codons is not a fixed difference.
        (["ATGTAA", "ATGTAA"], ["ATGTAG"], (0, 0, 0, 0)),
        (["ATGTAG", "ATGTAG"], ["ATGTGA"], (0, 0, 0, 0)),
        # A stop allele leaves only one clean allele, so the codon is not polymorphic.
        (["ATGTAC", "ATGTAA"], ["ATGTAC"], (0, 0, 0, 0)),
    ],
)
def test_stop_codons_add_no_differences(ingroup, outgroup, expected) -> None:
    assert counts(mk_test(sequence_set(ingroup), sequence_set(outgroup))) == expected


def test_terminal_stop_adds_no_sites() -> None:
    """Only ATG is counted. A stop codon has no synonymous sites, so it would add 3 to Ln."""
    result = mk_test(sequence_set(["ATGTAA", "ATGTAA"]), sequence_set(["ATGTAA"]))
    assert (result.ln, result.ls) == (pytest.approx(3.0), pytest.approx(0.0))


INTERNAL_INGROUP = ["ATGTAACCC", "ATGAAACCC", "ATGAAGCCC"]
INTERNAL_OUTGROUP = ["ATGAAACCC"]


def test_internal_stop_is_missing_data_for_that_sequence() -> None:
    result = mk_test(sequence_set(INTERNAL_INGROUP), sequence_set(INTERNAL_OUTGROUP))
    assert counts(result) == (0, 0, 0, 1)
    assert result.internal_stop_codons == 1


def test_internal_stop_leaves_the_frequency_denominator() -> None:
    """The stop sequence is not in the denominator: AAG is 1 of 2 clean copies."""
    data = extract_polymorphism_data(
        sequence_set(INTERNAL_INGROUP), sequence_set(INTERNAL_OUTGROUP)
    )
    assert data.polymorphisms == [(pytest.approx(0.5), "S")]
    assert data.internal_stop_codons == 1


def test_terminal_stop_is_not_an_internal_stop() -> None:
    result = mk_test(sequence_set(["ATGAAATAA", "ATGAAGTAA"]), sequence_set(["ATGAAATAA"]))
    assert result.internal_stop_codons == 0


class TestGeneticCode:
    """The genetic code of the analysis decides, not the standard code that the sets carry."""

    def test_code_stop_is_masked(self) -> None:
        """AGA is a stop codon under table 2, so AGA against CGA is not a difference."""
        result = mk_test(
            sequence_set(["ATGAGA", "ATGAGA"]),
            sequence_set(["ATGCGA"]),
            genetic_code=VERTEBRATE_MITO,
        )
        assert counts(result) == (0, 0, 0, 0)

    def test_standard_stop_is_data_under_another_code(self) -> None:
        """TGA is tryptophan under table 2, so TGA against TGG is a synonymous difference."""
        result = mk_test(
            sequence_set(["ATGTGA", "ATGTGA"]),
            sequence_set(["ATGTGG"]),
            genetic_code=VERTEBRATE_MITO,
        )
        assert counts(result) == (0, 1, 0, 0)


class TestPolarized:
    def test_stop_in_ingroup_is_not_a_fixed_difference(self) -> None:
        result = polarized_mk_test(
            sequence_set(["ATGTAA", "ATGTAA"]),
            sequence_set(["ATGTAG"]),
            sequence_set(["ATGTAG"]),
        )
        assert (result.dn_ingroup, result.ds_ingroup) == (0, 0)
        assert (result.dn_unpolarized, result.ds_unpolarized) == (0, 0)

    def test_stop_allele_is_not_a_polymorphism(self) -> None:
        result = polarized_mk_test(
            sequence_set(["ATGTAG", "ATGAAG"]),
            sequence_set(["ATGAAG"]),
            sequence_set(["ATGTAG"]),
        )
        assert (result.pn_ingroup, result.ps_ingroup) == (0, 0)

    def test_outgroup2_is_masked(self) -> None:
        pair = PolarizedAlignedPair(
            ingroup=sequence_set(["ATGAAA"]),
            outgroup=sequence_set(["ATGAAG"]),
            outgroup2=sequence_set(["ATGTAA"]),
        )
        assert pair.outgroup2.sequences[0].sequence == "ATGNNN"


class TestWithStopsMasked:
    def test_masks_each_stop_codon(self) -> None:
        original = sequence_set(["ATGTAACCC", "ATGAAATAA"])
        masked, stop_codons = original.with_stops_masked(DEFAULT_CODE)
        assert [s.sequence for s in masked.sequences] == ["ATGNNNCCC", "ATGAAANNN"]
        # Codon 2 is the terminal stop of the second sequence, so only codon 1 is internal.
        assert stop_codons == {1}
        assert masked.alignment_length == original.alignment_length

    def test_original_is_not_changed(self) -> None:
        original = sequence_set(["ATGTAACCC"])
        original.with_stops_masked(DEFAULT_CODE)
        assert original.sequences[0].sequence == "ATGTAACCC"

    def test_uses_the_reading_frame(self) -> None:
        """In frame 2 the codons are ATG, TAA and CCC, after the first base."""
        masked, stop_codons = sequence_set(["CATGTAACCC"], reading_frame=2).with_stops_masked(
            DEFAULT_CODE
        )
        assert masked.sequences[0].sequence == "CATGNNNCCC"
        assert masked.reading_frame == 2
        assert stop_codons == {1}


class TestInternalIsPerSequence:
    """A stop is terminal when its own sequence has no base after it, whatever the others have."""

    def test_stop_before_gap_padding_is_terminal(self) -> None:
        result = mk_test(
            sequence_set(["ATGAAATAA---", "ATGAAGTAA---"]), sequence_set(["ATGAAACTGCCC"])
        )
        assert result.internal_stop_codons == 0

    @pytest.mark.parametrize("padding", ["NNN", "???", "-N-"])
    def test_any_codon_with_no_base_is_padding(self, padding: str) -> None:
        _, stop_codons = sequence_set([f"ATGTAA{padding}"]).with_stops_masked(DEFAULT_CODE)
        assert stop_codons == set()

    def test_stop_with_a_base_later_is_internal(self) -> None:
        """Gaps after the stop do not end the sequence when a base comes later."""
        _, stop_codons = sequence_set(["ATGTAA---CCC"]).with_stops_masked(DEFAULT_CODE)
        assert stop_codons == {1}
