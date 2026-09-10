"""Builders for the hand-written alignments used across the test suite."""

from __future__ import annotations

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import Sequence, SequenceSet


def sequence_set(
    seqs: list[str], reading_frame: int = 1, genetic_code: GeneticCode = DEFAULT_CODE
) -> SequenceSet:
    """Build a SequenceSet from sequence strings, named s0, s1, and so on."""
    return SequenceSet(
        sequences=[Sequence(name=f"s{i}", sequence=s) for i, s in enumerate(seqs)],
        reading_frame=reading_frame,
        genetic_code=genetic_code,
    )
