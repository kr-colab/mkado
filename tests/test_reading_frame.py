"""The reading frame of SequenceSet arguments, and the reading_frame argument of each analysis."""

from __future__ import annotations

from collections.abc import Callable
from operator import attrgetter
from pathlib import Path

import pytest

from mkado.analysis.asymptotic import asymptotic_mk_test, extract_polymorphism_data
from mkado.analysis.mk_test import mk_test
from mkado.analysis.polarized import polarized_mk_test
from tests.builders import sequence_set, write_alignment_file

# In frame 2 the first base is skipped, so the codons are ATG AAA/AAG CCC.
INGROUP = ["CATGAAACCC", "CATGAAGCCC"]
OUTGROUP = ["CATGGAACCC"]
OUTGROUP2 = ["CATGAAACCC"]


def run_mk(frame: int, **kwargs):
    return mk_test(sequence_set(INGROUP, frame), sequence_set(OUTGROUP, frame), **kwargs)


def run_extract(frame: int, **kwargs):
    return extract_polymorphism_data(
        sequence_set(INGROUP, frame), sequence_set(OUTGROUP, frame), **kwargs
    )


def run_asymptotic(frame: int, **kwargs):
    return asymptotic_mk_test(
        sequence_set(INGROUP, frame),
        sequence_set(OUTGROUP, frame),
        bootstrap_replicates=0,
        **kwargs,
    )


def run_polarized(frame: int, **kwargs):
    return polarized_mk_test(
        sequence_set(INGROUP, frame),
        sequence_set(OUTGROUP, frame),
        sequence_set(OUTGROUP2, frame),
        **kwargs,
    )


ANALYSES: list[Callable] = [run_mk, run_extract, run_asymptotic, run_polarized]


@pytest.mark.parametrize("run", ANALYSES, ids=lambda f: f.__name__)
def test_explicit_mismatch_raises(run: Callable) -> None:
    """A reading_frame that differs from the frame of the sets would be ignored, so it is an error."""
    with pytest.raises(ValueError, match="frame"):
        run(1, reading_frame=2)


@pytest.mark.parametrize("run", ANALYSES, ids=lambda f: f.__name__)
def test_explicit_match_passes(run: Callable) -> None:
    run(2, reading_frame=2)


def test_sets_use_their_own_frame(tmp_path: Path) -> None:
    """Sets in frame 2 with no argument give the same counts as files read in frame 2."""
    ingroup = write_alignment_file(tmp_path / "in.fa", {f"s{i}": s for i, s in enumerate(INGROUP)})
    outgroup = write_alignment_file(tmp_path / "out.fa", {"s0": OUTGROUP[0]})
    from_files = mk_test(ingroup, outgroup, reading_frame=2)
    from_sets = run_mk(2)
    counts = attrgetter("dn", "ds", "pn", "ps")
    assert counts(from_sets) == counts(from_files)


def test_path_and_set_in_other_frame_raises(tmp_path: Path) -> None:
    """With no argument the path loads in frame 1, which does not match a set in frame 2."""
    ingroup = write_alignment_file(tmp_path / "in.fa", {f"s{i}": s for i, s in enumerate(INGROUP)})
    with pytest.raises(ValueError, match="frame"):
        mk_test(ingroup, sequence_set(OUTGROUP, 2))


def test_polarized_without_outgroup2_accepts_a_frame() -> None:
    """The second outgroup may be None; a reading_frame argument must not break that."""
    result = polarized_mk_test(
        sequence_set(INGROUP, reading_frame=1),
        sequence_set(OUTGROUP, reading_frame=1),
        None,
        reading_frame=1,
    )
    assert result is not None
