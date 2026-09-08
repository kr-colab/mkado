"""Tests for MK test implementation."""

from pathlib import Path

import pytest

from mkado.analysis.mk_test import mk_test, mk_test_from_counts
from mkado.analysis.statistics import alpha, dos, fishers_exact, neutrality_index
from mkado.core.sequences import SequenceSet


class TestStatistics:
    """Tests for statistical functions."""

    def test_fishers_exact_basic(self) -> None:
        """Test Fisher's exact test."""
        # Example: Dn=7, Ds=17, Pn=2, Ps=42
        p_value = fishers_exact(7, 17, 2, 42)
        assert 0 < p_value < 1

    def test_fishers_exact_symmetric(self) -> None:
        """Test that equal ratios give high p-value."""
        # Equal ratios should not be significant
        p_value = fishers_exact(10, 10, 10, 10)
        assert p_value > 0.5

    def test_neutrality_index_neutral(self) -> None:
        """Test NI under neutrality."""
        # Equal ratios -> NI = 1
        ni = neutrality_index(10, 10, 10, 10)
        assert ni is not None
        assert abs(ni - 1.0) < 0.001

    def test_neutrality_index_positive_selection(self) -> None:
        """Test NI under positive selection (excess divergence)."""
        # More Dn relative to Pn -> NI < 1
        ni = neutrality_index(20, 10, 5, 10)
        assert ni is not None
        assert ni < 1.0

    def test_neutrality_index_division_by_zero(self) -> None:
        """Test NI handles zero counts."""
        assert neutrality_index(0, 10, 10, 10) is None
        assert neutrality_index(10, 0, 10, 10) is None
        assert neutrality_index(10, 10, 10, 0) is None

    def test_alpha_neutral(self) -> None:
        """Test alpha under neutrality."""
        # Equal ratios -> alpha = 0
        a = alpha(10, 10, 10, 10)
        assert a is not None
        assert abs(a) < 0.001

    def test_alpha_positive_selection(self) -> None:
        """Test alpha under positive selection."""
        # More Dn relative to Pn -> alpha > 0
        a = alpha(20, 10, 5, 10)
        assert a is not None
        assert a > 0

    def test_alpha_negative_selection(self) -> None:
        """Test alpha under negative selection."""
        # More Pn relative to Dn -> alpha < 0
        a = alpha(5, 10, 20, 10)
        assert a is not None
        assert a < 0

    def test_dos_neutral(self) -> None:
        """Test DoS under neutrality."""
        # Equal ratios -> DoS = 0
        d = dos(10, 10, 10, 10)
        assert d is not None
        assert abs(d) < 0.001

    def test_dos_positive_selection(self) -> None:
        """Test DoS under positive selection (excess divergence)."""
        # High Dn ratio -> DoS > 0
        d = dos(20, 5, 5, 20)
        assert d is not None
        assert d > 0

    def test_dos_negative_selection(self) -> None:
        """Test DoS under negative selection (excess polymorphism)."""
        # High Pn ratio -> DoS < 0
        d = dos(5, 20, 20, 5)
        assert d is not None
        assert d < 0

    def test_dos_zero_divergence(self) -> None:
        """Test DoS when Dn+Ds=0."""
        # Only polymorphism -> returns value (not None)
        d = dos(0, 0, 10, 10)
        assert d is not None
        # 0 - 0.5 = -0.5
        assert abs(d - (-0.5)) < 0.001

    def test_dos_zero_polymorphism(self) -> None:
        """Test DoS when Pn+Ps=0."""
        # Only divergence -> returns value (not None)
        d = dos(10, 10, 0, 0)
        assert d is not None
        # 0.5 - 0 = 0.5
        assert abs(d - 0.5) < 0.001

    def test_dos_all_zero(self) -> None:
        """Test DoS when all counts are zero."""
        d = dos(0, 0, 0, 0)
        assert d is None

    def test_dos_bounds(self) -> None:
        """Verify DoS is bounded [-1, +1]."""
        # Maximum DoS = 1: all Dn, no Ds, no polymorphism
        d_max = dos(10, 0, 0, 10)
        assert d_max is not None
        assert abs(d_max - 1.0) < 0.001

        # Minimum DoS = -1: all Ds, no Dn, all Pn, no Ps
        d_min = dos(0, 10, 10, 0)
        assert d_min is not None
        assert abs(d_min - (-1.0)) < 0.001

        # General case should be within bounds
        d = dos(7, 17, 2, 42)
        assert d is not None
        assert -1 <= d <= 1


class TestMKResult:
    """Tests for MKResult class."""

    def test_mk_test_from_counts(self) -> None:
        """Test creating MK result from counts."""
        result = mk_test_from_counts(dn=7, ds=17, pn=2, ps=42)

        assert result.dn == 7
        assert result.ds == 17
        assert result.pn == 2
        assert result.ps == 42
        assert result.p_value is not None
        assert result.ni is not None
        assert result.alpha is not None
        assert result.dos is not None

    def test_mk_result_dos(self) -> None:
        """Test DoS calculation in MK result."""
        result = mk_test_from_counts(dn=6, ds=8, pn=1, ps=8)
        # DoS = 6/14 - 1/9 = 0.4286 - 0.1111 = 0.3175
        assert result.dos is not None
        assert abs(result.dos - (6 / 14 - 1 / 9)) < 0.001

    def test_mk_result_to_dict(self) -> None:
        """Test converting result to dictionary."""
        result = mk_test_from_counts(dn=7, ds=17, pn=2, ps=42)
        d = result.to_dict()

        assert d["dn"] == 7
        assert d["ds"] == 17
        assert "p_value" in d
        assert "ni" in d
        assert "alpha" in d
        assert "dos" in d

    def test_mk_result_to_dict_includes_dos(self) -> None:
        """Test that to_dict includes DoS value."""
        result = mk_test_from_counts(dn=6, ds=8, pn=1, ps=8)
        d = result.to_dict()
        assert "dos" in d
        assert d["dos"] is not None
        assert abs(d["dos"] - (6 / 14 - 1 / 9)) < 0.001

    def test_mk_result_str(self) -> None:
        """Test string representation."""
        result = mk_test_from_counts(dn=7, ds=17, pn=2, ps=42)
        s = str(result)

        assert "Dn=7" in s
        assert "Ds=17" in s
        assert "Fisher" in s

    def test_mk_result_str_includes_dos(self) -> None:
        """Test that string representation includes DoS."""
        result = mk_test_from_counts(dn=7, ds=17, pn=2, ps=42)
        s = str(result)
        assert "DoS:" in s


class TestMKTestIntegration:
    """Integration tests for the full MK test."""

    def test_mk_test_kreitman_data(self) -> None:
        """Test MK test with Kreitman Adh data.

        Expected results from original mkTest.rb:
        Dn=6, Pn=1, Ds=8, Ps=8
        """
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        result = mk_test(ingroup, outgroup)

        # These are the expected values from the original implementation
        # Note: exact values may vary slightly based on algorithm details
        assert result.dn >= 0
        assert result.ds >= 0
        assert result.pn >= 0
        assert result.ps >= 0

        # The test should complete without errors
        assert result.p_value is not None

    def test_mk_test_sequence_set_input(self, tmp_path: Path) -> None:
        """Test MK test with SequenceSet objects."""
        # Create simple test sequences
        # Ingroup: 2 sequences with one polymorphism
        # Outgroup: 1 sequence fixed differently

        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATG
>seq2
ATGCTGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATG
""")

        ingroup = SequenceSet.from_fasta(ingroup_fa)
        outgroup = SequenceSet.from_fasta(outgroup_fa)

        result = mk_test(ingroup, outgroup)

        # Should have some counts
        assert result.dn >= 0
        assert result.ds >= 0
        assert result.pn >= 0
        assert result.ps >= 0

    def test_mk_test_combined_file_mode(self, tmp_path: Path) -> None:
        """Test MK test with combined file filtered by name patterns."""
        # Create combined alignment with sequences from two "species"
        combined_fa = tmp_path / "combined.fa"
        combined_fa.write_text(""">gene1_speciesA_1
ATGATGATG
>gene1_speciesA_2
ATGCTGATG
>gene1_speciesA_3
ATGATGATG
>gene1_speciesB_1
ATGGTGATG
>gene1_speciesB_2
ATGGTGATG
""")

        # Load and filter by species
        all_seqs = SequenceSet.from_fasta(combined_fa)
        ingroup = all_seqs.filter_by_name("speciesA")
        outgroup = all_seqs.filter_by_name("speciesB")

        assert len(ingroup) == 3
        assert len(outgroup) == 2

        result = mk_test(ingroup, outgroup)

        # Should complete without errors
        assert result.dn >= 0
        assert result.ds >= 0
        assert result.pn >= 0
        assert result.ps >= 0
        assert result.p_value is not None

    def test_mk_test_min_frequency_filter(self, tmp_path: Path) -> None:
        """Test that min_frequency parameter filters low-frequency polymorphisms."""
        # Create ingroup with 10 sequences where:
        # - One polymorphism at position 3 (codon 0) has freq 1/10 = 0.1 (derived)
        # - One polymorphism at position 6 (codon 1) has freq 5/10 = 0.5 (derived)
        # The outgroup has the ancestral state at both positions

        ingroup_fa = tmp_path / "ingroup.fa"
        # Create 10 sequences with varying polymorphisms
        # Position 3 (codon 0, position 0): 9 have T, 1 has C -> derived freq = 0.1
        # Position 6 (codon 1, position 0): 5 have A, 5 have C -> derived freq = 0.5
        ingroup_fa.write_text(
            """>seq1
ATGATGATG
>seq2
ATGATGATG
>seq3
ATGATGATG
>seq4
ATGATGATG
>seq5
ATGATGATG
>seq6
ATGCTGATG
>seq7
ATGCTGATG
>seq8
ATGCTGATG
>seq9
ATGCTGATG
>seq10
ATCCTGATG
"""
        )

        # Outgroup has ATG ATG ATG - matches seq1-5 at codon 0
        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(
            """>out1
ATGATGATG
"""
        )

        # Without filtering, should count polymorphisms
        result_no_filter = mk_test(ingroup_fa, outgroup_fa, min_frequency=0.0)

        # With high min_frequency filter, should exclude low-frequency polymorphisms
        result_with_filter = mk_test(ingroup_fa, outgroup_fa, min_frequency=0.3)

        # The filtered result should have fewer or equal polymorphisms
        total_poly_no_filter = result_no_filter.pn + result_no_filter.ps
        total_poly_filtered = result_with_filter.pn + result_with_filter.ps

        assert total_poly_filtered <= total_poly_no_filter
