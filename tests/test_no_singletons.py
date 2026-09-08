"""--no-singletons on the FASTA analysis paths.

A singleton is one copy of the derived allele, whatever the sample size, whether some
sequences are missing data at that codon, and whether polymorphisms are pooled. Each
path is checked against all three, because a frequency threshold cannot express them.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mkado.analysis.asymptotic import extract_polymorphism_data
from mkado.analysis.mk_test import mk_test
from mkado.analysis.polarized import polarized_mk_test


# Codon 1 of ATG GCC AAA carries the variation. GCC and GCT are both Ala, so a
# difference there is synonymous, and the outgroup holds the ancestral GCC.
ANCESTRAL = "ATGGCCAAA"
DERIVED = "ATGGCTAAA"


def singleton_alignment(n: int = 5) -> dict[str, str]:
    """One ingroup sequence carries the derived codon; the rest are ancestral."""
    records = {f"ingroup_{i}": ANCESTRAL for i in range(1, n)}
    records[f"ingroup_{n}"] = DERIVED
    return records


@pytest.fixture
def singleton(tmp_path: Path, write_alignment):
    """Five ingroup sequences, one carrying a synonymous derived codon."""
    ingroup = write_alignment(tmp_path / "in.fa", singleton_alignment(5))
    outgroup = write_alignment(tmp_path / "out.fa", {"outgroup_1": ANCESTRAL})
    return ingroup, outgroup


def _standard(ingroup, outgroup, **kwargs) -> tuple[int, int]:
    result = mk_test(ingroup, outgroup, **kwargs)
    return result.pn, result.ps


def _polarized(ingroup, outgroup, outgroup2, **kwargs) -> tuple[int, int]:
    result = polarized_mk_test(ingroup, outgroup, outgroup2, **kwargs)
    return result.pn_ingroup, result.ps_ingroup


def _extracted(ingroup, outgroup, **kwargs) -> tuple[int, int]:
    data = extract_polymorphism_data(ingroup, outgroup, **kwargs)
    nonsyn = sum(1 for _, kind in data.polymorphisms if kind == "N")
    return nonsyn, len(data.polymorphisms) - nonsyn


class TestStandard:
    def test_singleton_is_kept_by_default(self, singleton):
        assert _standard(*singleton) == (0, 1)

    @pytest.mark.parametrize(
        ("derived", "expected"),
        [("ATGGCTAAA", (0, 1)), ("ATGCCCAAA", (1, 0))],
        ids=["synonymous", "nonsynonymous"],
    )
    def test_singleton_is_removed(self, tmp_path, write_alignment, derived, expected):
        records = singleton_alignment(5)
        records["ingroup_5"] = derived
        ingroup = write_alignment(tmp_path / "in.fa", records)
        outgroup = write_alignment(tmp_path / "out.fa", {"outgroup_1": ANCESTRAL})
        assert _standard(ingroup, outgroup) == expected
        assert _standard(ingroup, outgroup, no_singletons=True) == (0, 0)

    def test_two_carriers_are_not_a_singleton(self, tmp_path, write_alignment):
        records = singleton_alignment(5)
        records["ingroup_4"] = DERIVED
        ingroup = write_alignment(tmp_path / "in.fa", records)
        outgroup = write_alignment(tmp_path / "out.fa", {"outgroup_1": ANCESTRAL})
        assert _standard(ingroup, outgroup, no_singletons=True) == (0, 1)

    def test_singleton_removed_when_another_sequence_has_a_gap(self, tmp_path, write_alignment):
        """The gap shrinks the denominator, so the frequency is no longer one over n."""
        records = singleton_alignment(5)
        records["ingroup_1"] = "ATG---AAA"
        ingroup = write_alignment(tmp_path / "in.fa", records)
        outgroup = write_alignment(tmp_path / "out.fa", {"outgroup_1": ANCESTRAL})
        assert _standard(ingroup, outgroup) == (0, 1)
        assert _standard(ingroup, outgroup, no_singletons=True) == (0, 0)

    def test_singleton_removed_when_pooling(self, singleton):
        """Pooling adds outgroup sequences, which must not change what a singleton is."""
        ingroup, outgroup = singleton
        assert _standard(ingroup, outgroup, pool_polymorphisms=True, no_singletons=True) == (0, 0)

    def test_pooled_outgroup_only_site_is_not_a_singleton(self, tmp_path, write_alignment):
        """The ingroup is invariant here, so it contributes no derived copies at all.

        Counting zero copies is not counting one, and the outgroup carries three.
        """
        ingroup = write_alignment(tmp_path / "in.fa", {f"i{i}": ANCESTRAL for i in range(1, 4)})
        outgroup = write_alignment(
            tmp_path / "out.fa",
            {"o1": DERIVED, "o2": ANCESTRAL, "o3": ANCESTRAL, "o4": ANCESTRAL},
        )
        kept = _standard(ingroup, outgroup, pool_polymorphisms=True)
        assert _standard(ingroup, outgroup, pool_polymorphisms=True, no_singletons=True) == kept

    def test_site_at_min_frequency_is_kept(self, singleton):
        """Only sites below the minimum are dropped, as the option's documentation says."""
        assert _standard(*singleton, min_frequency=0.2) == (0, 1)
        assert _standard(*singleton, min_frequency=0.3) == (0, 0)


class TestPolarized:
    def test_singleton_is_kept_by_default(self, singleton):
        ingroup, outgroup = singleton
        assert _polarized(ingroup, outgroup, outgroup) == (0, 1)

    def test_singleton_is_removed(self, singleton):
        ingroup, outgroup = singleton
        assert _polarized(ingroup, outgroup, outgroup, no_singletons=True) == (0, 0)

    def test_site_at_min_frequency_is_kept(self, singleton):
        ingroup, outgroup = singleton
        assert _polarized(ingroup, outgroup, outgroup, min_frequency=0.2) == (0, 1)


class TestExtractedForAlphaTg:
    def test_singleton_is_kept_by_default(self, singleton):
        assert _extracted(*singleton) == (0, 1)

    def test_singleton_is_removed(self, singleton):
        assert _extracted(*singleton, no_singletons=True) == (0, 0)

    def test_site_at_min_frequency_is_kept(self, singleton):
        assert _extracted(*singleton, min_frequency=0.2) == (0, 1)
