"""Tests for Sequence and SequenceSet."""

import pickle

import numpy as np
import pytest

from mkado.core.codons import GeneticCode
from mkado.core.sequences import AlignmentLengthError, Sequence, SequenceSet
from tests.builders import sequence_set, write_alignment_file

VERTEBRATE_MITO = GeneticCode(table_id=2)

# Three codons per sequence. Codon 1 is polymorphic and one sequence has a gap there.
ALIGNMENT = ["ATGAAACCC", "ATGAAGCCC", "ATG---CCC"]


class TestSequence:
    """Accessors on a single sequence."""

    def test_len(self) -> None:
        """Length is the number of nucleotides."""
        assert len(Sequence("x", "ATGAAACCC")) == 9

    def test_getitem_index_and_slice(self) -> None:
        """Indexing and slicing read the underlying string."""
        seq = Sequence("x", "ATGAAACCC")
        assert seq[0] == "A"
        assert seq[3:6] == "AAA"

    @pytest.mark.parametrize(
        "reading_frame, codon_index, expected",
        [(1, 1, "AAA"), (2, 0, "TGA"), (3, 0, "GAA")],
    )
    def test_get_codon_by_reading_frame(
        self, reading_frame: int, codon_index: int, expected: str
    ) -> None:
        """The reading frame offsets where codons start."""
        assert Sequence("x", "ATGAAACCC").get_codon(codon_index, reading_frame) == expected

    @pytest.mark.parametrize(
        "sequence, reading_frame, expected",
        [("ATGAAACCC", 1, 3), ("ATGAAACCC", 2, 2), ("ATGAAACCC", 3, 2), ("GATGAAACCC", 2, 3)],
    )
    def test_num_codons_by_reading_frame(
        self, sequence: str, reading_frame: int, expected: int
    ) -> None:
        """Only complete codons after the frame offset are counted."""
        assert Sequence("x", sequence).num_codons(reading_frame) == expected


class TestSequenceSetBasics:
    """Length, indexing, and the two size properties."""

    def test_len_and_getitem(self) -> None:
        """Indexing returns the Sequence at that position."""
        s = sequence_set(ALIGNMENT)
        assert len(s) == 3
        assert s[1].name == "s1"
        assert s[1].sequence == "ATGAAGCCC"

    def test_empty_set(self) -> None:
        """An empty set reports zero sizes instead of failing on the first sequence."""
        empty = SequenceSet()
        assert empty.num_codons == 0
        assert empty.alignment_length == 0
        assert empty.codon_array().shape == (0, 0)

    def test_alignment_length_is_first_sequence_length(self) -> None:
        """The alignment length comes from the first sequence."""
        assert sequence_set(ALIGNMENT).alignment_length == 9

    def test_num_codons_follows_reading_frame(self) -> None:
        """A leading base outside the frame does not count toward a codon."""
        assert sequence_set(["GATGAAACCC"], reading_frame=2).num_codons == 3


class TestPickling:
    """The codon cache is dropped on pickling and rebuilt on demand."""

    def test_pickle_drops_codon_cache(self) -> None:
        """A populated cache adds nothing to the pickle."""
        populated = sequence_set(ALIGNMENT)
        populated.codon_set_clean(0)
        populated.codon_set_clean(1)
        fresh = sequence_set(ALIGNMENT)

        assert len(pickle.dumps(populated)) == len(pickle.dumps(fresh))

    def test_round_trip_preserves_fields_and_answers(self) -> None:
        """The restored set has the same fields, gives the same codon sets, and caches again."""
        original = sequence_set(ALIGNMENT)
        original.codon_set_clean(1)

        restored = pickle.loads(pickle.dumps(original))

        # GeneticCode compares by identity, so the whole set never compares equal
        # after a round trip. Compare the fields instead.
        assert restored.sequences == original.sequences
        assert restored.reading_frame == original.reading_frame
        assert restored.genetic_code.code == original.genetic_code.code
        first = restored.codon_set_clean(1)
        assert first == original.codon_set_clean(1) == {"AAA", "AAG"}
        assert restored.codon_set_clean(1) is first

    def test_round_trip_keeps_non_default_genetic_code(self) -> None:
        """A non-standard genetic code survives the round trip."""
        original = sequence_set(["ATGAGA"], genetic_code=VERTEBRATE_MITO)

        restored = pickle.loads(pickle.dumps(original))

        assert restored.genetic_code.code == VERTEBRATE_MITO.code
        assert restored.amino_set(1) == {"*"}


class TestCodonAndSiteSets:
    """Per-position codon, amino acid, and nucleotide sets."""

    def test_codon_set_includes_missing_data(self) -> None:
        """codon_set keeps gapped codons and codon_set_clean drops them."""
        s = sequence_set(ALIGNMENT)
        assert s.codon_set(1) == {"AAA", "AAG", "---"}
        assert s.codon_set_clean(1) == {"AAA", "AAG"}

    def test_amino_set_translates_every_codon(self) -> None:
        """A gapped codon translates to X and is absent from the clean set."""
        s = sequence_set(ALIGNMENT)
        assert s.amino_set(1) == {"K", "X"}
        assert s.amino_set_clean(1) == {"K"}

    def test_amino_set_uses_the_genetic_code(self) -> None:
        """AGA is arginine under the standard code and a stop under table 2."""
        standard = sequence_set(["AGA"])
        mito = sequence_set(["AGA"], genetic_code=VERTEBRATE_MITO)
        assert standard.amino_set(0) == {"R"}
        assert mito.amino_set(0) == {"*"}

    def test_site_set_and_clean(self) -> None:
        """site_set keeps gaps, site_set_clean drops them, and a site past the end is empty."""
        s = sequence_set(ALIGNMENT)
        assert s.site_set(3) == {"A", "-"}
        assert s.site_set_clean(3) == {"A"}
        assert s.site_set(99) == set()


class TestCodonArray:
    """The sequences-by-codons array."""

    def test_shape_dtype_and_values(self) -> None:
        """One row per sequence, one column per codon, three-character strings."""
        arr = sequence_set(ALIGNMENT).codon_array()
        assert arr.shape == (3, 3)
        assert arr.dtype == np.dtype("U3")
        assert arr.tolist() == [
            ["ATG", "AAA", "CCC"],
            ["ATG", "AAG", "CCC"],
            ["ATG", "---", "CCC"],
        ]

    def test_reading_frame_offsets_codons(self) -> None:
        """Codons are read from the frame offset."""
        arr = sequence_set(["GATGAAACCC"], reading_frame=2).codon_array()
        assert arr.tolist() == [["ATG", "AAA", "CCC"]]


class TestFrequencies:
    """Codon counts and frequencies at a position."""

    def test_site_codon_counts_skips_missing_data(self) -> None:
        """A gapped sequence is left out of both the counts and the total."""
        assert sequence_set(ALIGNMENT).site_codon_counts(1) == ({"AAA": 1, "AAG": 1}, 2)

    def test_site_frequency_spectrum(self) -> None:
        """Frequencies divide by the sequences with data at the codon."""
        assert sequence_set(ALIGNMENT).site_frequency_spectrum(1) == {"AAA": 0.5, "AAG": 0.5}

    def test_site_frequency_spectrum_without_data_is_empty(self) -> None:
        """An all-gap column has no frequencies."""
        assert sequence_set(["---", "---"]).site_frequency_spectrum(0) == {}

    def test_polymorphic_codons(self) -> None:
        """Only the codon with two clean alleles is polymorphic."""
        assert sequence_set(ALIGNMENT).polymorphic_codons() == [1]


class TestDerivedState:
    """Derived allele frequency and count against an ancestral source."""

    OUTGROUP = ["ATGAAACCC"]

    def derived(
        self, ingroup: list[str], outgroup: list[str] | None = None, codon_index: int = 1
    ) -> tuple[float, int] | None:
        return sequence_set(ingroup).derived_state(
            sequence_set(outgroup or self.OUTGROUP), codon_index
        )

    def test_ancestral_present(self) -> None:
        """Every non-ancestral codon counts toward the derived total."""
        ingroup = ["ATGAAACCC", "ATGAAACCC", "ATGAAGCCC", "ATGAAGCCC", "ATGAATCCC"]
        assert self.derived(ingroup) == (0.6, 3)

    def test_invariant_site(self) -> None:
        """A site fixed for the ancestral codon has no derived copies."""
        assert self.derived(["ATGAAACCC", "ATGAAACCC"]) == (0.0, 0)

    def test_ancestral_absent(self) -> None:
        """No shared codon means the ancestral state is unknown."""
        assert self.derived(["ATGAAGCCC", "ATGAATCCC"]) is None

    def test_all_ingroup_missing(self) -> None:
        """A codon with no clean ingroup data has no derived state."""
        assert self.derived(["ATG---CCC", "ATGNNNCCC"]) is None

    def test_source_without_clean_codon(self) -> None:
        """An ancestral source with only missing data cannot polarize."""
        assert self.derived(["ATGAAACCC", "ATGAAGCCC"], outgroup=["ATG---CCC"]) is None

    def test_ancestral_is_most_frequent_shared_codon(self) -> None:
        """The ancestral codon is the shared one, even when another codon is commoner."""
        ingroup = ["ATGAAACCC", "ATGAAACCC", "ATGAAACCC", "ATGAAGCCC", "ATGAAGCCC"]
        assert self.derived(ingroup, outgroup=["ATGAAGCCC"]) == (0.6, 3)

    def test_missing_data_shrinks_denominator(self) -> None:
        """A gapped sequence is left out of the frequency denominator."""
        result = self.derived(["ATGAAACCC", "ATGAAGCCC", "ATG---CCC", "ATGAAACCC"])
        assert result is not None
        assert result[0] == pytest.approx(1 / 3)
        assert result[1] == 1

    def test_tie_between_shared_codons_is_deterministic(self) -> None:
        """Two shared codons at equal frequency give the same answer whichever is chosen."""
        ingroup = ["ATGAAACCC", "ATGAAACCC", "ATGAAGCCC", "ATGAAGCCC"]
        assert self.derived(ingroup, outgroup=["ATGAAACCC", "ATGAAGCCC"]) == (0.5, 2)

    def test_codon_index_selects_the_site(self) -> None:
        """The answer follows the codon index, not the first polymorphic site."""
        ingroup = ["ATGAAACCC", "ATGAAACCG"]
        assert self.derived(ingroup, codon_index=1) == (0.0, 0)
        assert self.derived(ingroup, codon_index=2) == (0.5, 1)


class TestAmbiguousBases:
    """A codon with any base other than A, C, G or T is missing data for that sequence."""

    @pytest.mark.parametrize("code", list("RYSWKMBDHV?."))
    def test_codon_set_clean_drops_ambiguous_codon(self, code: str) -> None:
        """An ambiguity code removes the codon, as N and gaps do."""
        s = sequence_set(["ATGAAA", f"ATGA{code}A", "ATGAAG"])
        assert s.codon_set_clean(1) == {"AAA", "AAG"}

    def test_site_codon_counts_drops_ambiguous_codon(self) -> None:
        """An ambiguous codon is left out of both the counts and the total."""
        s = sequence_set(["ATGAAA", "ATGARA", "ATGAAG", "ATGAYG"])
        assert s.site_codon_counts(1) == ({"AAA": 1, "AAG": 1}, 2)

    def test_site_set_clean_drops_ambiguous_base(self) -> None:
        """Only A, C, G and T remain in the clean nucleotide set."""
        s = sequence_set(["ATGAAA", "ATGRAA", "ATG?AA", "ATGCAA"])
        assert s.site_set_clean(3) == {"A", "C"}

    def test_lowercase_input_is_uppercased(self) -> None:
        """Lowercase bases are read as their uppercase form, not as missing data."""
        assert Sequence("x", "atgaaa").sequence == "ATGAAA"
        s = sequence_set(["atgaaa", "ATGAAG"])
        assert s.codon_set_clean(1) == {"AAA", "AAG"}

    def test_derived_state_denominator_drops_ambiguous_codon(self) -> None:
        """A sequence with an ambiguous codon is left out of the frequency denominator."""
        ingroup = sequence_set(["ATGAAACCC", "ATGAAGCCC", "ATGARACCC", "ATGAAACCC"])
        outgroup = sequence_set(["ATGAAACCC"])
        assert ingroup.derived_state(outgroup, 1) == (pytest.approx(1 / 3), 1)


class TestPositionDerivedState:
    """Polarization reads one codon position, so a difference elsewhere does not block it."""

    INGROUP = ["AAA", "AAG", "AAA", "AAG"]

    def state(self, ingroup: list[str], outgroup: list[str], position: int = 2):
        return sequence_set(ingroup).position_derived_state(sequence_set(outgroup), 0, position)

    def test_outgroup_differs_at_another_position(self) -> None:
        """CAA shares the A at position 2 with half the ingroup."""
        assert self.state(self.INGROUP, ["CAA"]) == (0.5, 2)

    def test_third_allele_at_the_position(self) -> None:
        assert self.state(self.INGROUP, ["AAC"]) is None

    def test_outgroup_with_no_clean_codon(self) -> None:
        assert self.state(self.INGROUP, ["NNN"]) is None

    @pytest.mark.parametrize("outgroup", [["AAA", "AAG"], ["AAG", "AAA"]])
    def test_tie_goes_to_the_first_base(self, outgroup: list[str]) -> None:
        """A and G are both shared and equally common: A is ancestral, whatever the order."""
        assert self.state(self.INGROUP, outgroup) == (0.5, 2)

    def test_missing_ingroup_codon_is_not_counted(self) -> None:
        assert self.state(["AAA", "AAG", "NNN", "AAA"], ["AAA"]) == (pytest.approx(1 / 3), 1)


class TestLengthCheck:
    """Every sequence in a set must have the same length."""

    def test_equal_lengths_load(self) -> None:
        assert sequence_set(ALIGNMENT).alignment_length == 9

    def test_unequal_lengths_raise(self) -> None:
        with pytest.raises(AlignmentLengthError) as excinfo:
            sequence_set(["ATGAAA", "ATG"])
        assert isinstance(excinfo.value, ValueError)

    def test_error_does_not_depend_on_order(self) -> None:
        """The two orderings from issue 136 give the same lengths in the message."""
        with pytest.raises(AlignmentLengthError) as short_last:
            sequence_set(["ATGAAA", "ATGAAG", "ATG"])
        with pytest.raises(AlignmentLengthError) as short_first:
            sequence_set(["ATG", "ATGAAG", "ATGAAA"])
        for excinfo in (short_last, short_first):
            assert "3 bp" in str(excinfo.value)
            assert "6 bp" in str(excinfo.value)

    def test_one_base_insertion_raises(self) -> None:
        with pytest.raises(AlignmentLengthError):
            sequence_set(["ATGAAA", "ATGAAAA"])

    def test_from_fasta_names_file_and_lengths(self, tmp_path) -> None:
        fasta = write_alignment_file(
            tmp_path / "ragged.fa", {"in3": "ATG", "in2": "ATGAAG", "in1": "ATGAAA"}
        )
        with pytest.raises(AlignmentLengthError) as excinfo:
            SequenceSet.from_fasta(fasta)
        message = str(excinfo.value)
        assert "ragged.fa" in message
        assert "3 bp" in message
        assert "6 bp" in message
        assert "in3" in message

    def test_long_name_list_is_cut(self) -> None:
        """A length shared by many sequences names only the first few."""
        with pytest.raises(AlignmentLengthError) as excinfo:
            sequence_set(["ATGAAA"] + ["ATG"] * 10)
        message = str(excinfo.value)
        assert "and 7 more" in message
        assert "s10" not in message
