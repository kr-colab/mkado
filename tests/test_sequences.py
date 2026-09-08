"""Focused tests for SequenceSet array, frequency, and serialization helpers."""

from __future__ import annotations

import pickle

import numpy as np
import pytest

from mkado.core.sequences import Sequence, SequenceSet


def sequence_set(*sequences: str) -> SequenceSet:
    return SequenceSet(
        sequences=[Sequence(name=f"seq{index}", sequence=value) for index, value in enumerate(sequences)]
    )


def test_pickle_round_trip_drops_and_rebuilds_codon_cache() -> None:
    original = sequence_set("AAATTT", "AAGTTT")
    expected = original.codon_set_clean(0)
    assert original._codon_set_clean_cache == {0: expected}

    restored = pickle.loads(pickle.dumps(original))

    assert restored.sequences == original.sequences
    assert restored.reading_frame == original.reading_frame
    assert restored.genetic_code.code == original.genetic_code.code
    assert restored._codon_set_clean_cache == {}
    assert restored.codon_set_clean(0) == expected
    assert restored._codon_set_clean_cache == {0: expected}


def test_empty_sequence_set_has_zero_alignment_length_and_empty_codon_array() -> None:
    sequences = SequenceSet()
    assert sequences.alignment_length == 0
    assert sequences.codon_array().shape == (0, 0)


def test_codon_array_has_expected_shape_dtype_and_values() -> None:
    sequences = sequence_set(
        "AAATTTCCCGGG",
        "AAGTTTCCCGGA",
        "AAATTCCTCGGG",
    )

    result = sequences.codon_array()

    assert result.shape == (3, 4)
    assert result.dtype == np.dtype("<U3")
    np.testing.assert_array_equal(
        result,
        np.array(
            [
                ["AAA", "TTT", "CCC", "GGG"],
                ["AAG", "TTT", "CCC", "GGA"],
                ["AAA", "TTC", "CTC", "GGG"],
            ]
        ),
    )


@pytest.mark.parametrize(
    ("sequences", "ancestral_codon", "expected"),
    [
        (("AAA", "AAA", "AAG"), "AAA", pytest.approx(1 / 3)),
        (("AAA", "AAG", "AAG"), "CCC", 1.0),
        (("AAA", "AAA", "AAA"), "aaa", 0.0),
        (("NNN", "---"), "AAA", None),
    ],
)
def test_derived_allele_frequency(
    sequences: tuple[str, ...], ancestral_codon: str, expected: float | None
) -> None:
    assert sequence_set(*sequences).derived_allele_frequency(0, ancestral_codon) == expected
