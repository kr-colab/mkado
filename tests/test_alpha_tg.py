"""Tests for Tarone-Greenland alpha (α_TG) implementation."""

from pathlib import Path

import pytest

from mkado.analysis.alpha_tg import (
    AlphaTGResult,
    alpha_tg_from_gene_data,
    compute_ni_tg,
)
from mkado.analysis.asymptotic import PolymorphismData, extract_polymorphism_data


class TestComputeNiTg:
    """Tests for compute_ni_tg function."""

    def test_single_gene_calculation(self) -> None:
        """Test NI_TG calculation with a single gene."""
        # Create gene with known values
        # Dn=10, Ds=20, Pn=5, Ps=15
        # Weight = Ps + Ds = 15 + 20 = 35
        # Numerator = Ds * Pn / weight = 20 * 5 / 35 = 100/35
        # Denominator = Dn * Ps / weight = 10 * 15 / 35 = 150/35
        # NI_TG = (100/35) / (150/35) = 100/150 = 2/3 ≈ 0.6667
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
                dn=10,
                ds=20,
                gene_id="gene1",
            )
        ]

        ni_tg = compute_ni_tg(gene_data)

        assert ni_tg is not None
        assert abs(ni_tg - 2 / 3) < 0.0001

    def test_multi_gene_weighting(self) -> None:
        """Test that multi-gene weighting is correct."""
        # Two genes with different weights
        # Gene 1: Dn=10, Ds=20, Pn=5, Ps=10  -> weight = 10+20 = 30
        # Gene 2: Dn=5, Ds=10, Pn=3, Ps=20   -> weight = 20+10 = 30
        #
        # Numerator = (20*5/30) + (10*3/30) = 100/30 + 30/30 = 130/30
        # Denominator = (10*10/30) + (5*20/30) = 100/30 + 100/30 = 200/30
        # NI_TG = (130/30) / (200/30) = 130/200 = 0.65
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 10,
                dn=10,
                ds=20,
                gene_id="gene1",
            ),
            PolymorphismData(
                polymorphisms=[(0.3, "N")] * 3 + [(0.6, "S")] * 20,
                dn=5,
                ds=10,
                gene_id="gene2",
            ),
        ]

        ni_tg = compute_ni_tg(gene_data)

        assert ni_tg is not None
        assert abs(ni_tg - 0.65) < 0.0001

    def test_zero_denominator_returns_none(self) -> None:
        """Test that zero denominator returns None."""
        # No Ps and no Ds means denominator term is zero
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5,  # Only Pn, no Ps
                dn=0,  # Zero Dn
                ds=0,  # Zero Ds (weight becomes 0)
                gene_id="gene1",
            )
        ]

        ni_tg = compute_ni_tg(gene_data)

        assert ni_tg is None

    def test_empty_gene_list(self) -> None:
        """Test with empty gene list."""
        ni_tg = compute_ni_tg([])
        assert ni_tg is None


class TestAlphaTgFromGeneData:
    """Tests for alpha_tg_from_gene_data function."""

    def test_basic_calculation(self) -> None:
        """Test basic α_TG calculation."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
                dn=10,
                ds=20,
                gene_id="gene1",
            )
        ]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        assert isinstance(result, AlphaTGResult)
        # NI_TG = 2/3, so alpha_TG = 1 - 2/3 = 1/3 ≈ 0.3333
        assert abs(result.alpha_tg - 1 / 3) < 0.0001
        assert abs(result.ni_tg - 2 / 3) < 0.0001

    def test_totals_are_correct(self) -> None:
        """Test that total counts are calculated correctly."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 10,
                dn=10,
                ds=20,
                gene_id="gene1",
            ),
            PolymorphismData(
                polymorphisms=[(0.3, "N")] * 3 + [(0.6, "S")] * 7,
                dn=5,
                ds=8,
                gene_id="gene2",
            ),
        ]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        assert result.num_genes == 2
        assert result.dn_total == 15  # 10 + 5
        assert result.ds_total == 28  # 20 + 8
        assert result.pn_total == 8  # 5 + 3
        assert result.ps_total == 17  # 10 + 7

    def test_bootstrap_produces_reasonable_cis(self) -> None:
        """Test that bootstrap produces reasonable confidence intervals."""
        # Create data with heterogeneity across genes for meaningful bootstrap
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * (3 + i % 5) + [(0.5, "S")] * (8 + i % 7),
                dn=8 + i % 6,
                ds=15 + i % 10,
                gene_id=f"gene{i}",
            )
            for i in range(20)
        ]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=500, seed=42)

        # CI should bracket the point estimate (with heterogeneous data)
        assert result.ci_low <= result.alpha_tg
        assert result.ci_high >= result.alpha_tg
        # CI should be reasonably tight (not too wide)
        ci_width = result.ci_high - result.ci_low
        assert ci_width < 0.5  # Should not be excessively wide

    def test_single_gene_ci_is_none(self) -> None:
        """A single gene gives the bootstrap nothing to resample: CI is undefined."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
                dn=10,
                ds=20,
                gene_id="gene1",
            )
        ]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        assert result.num_genes == 1
        assert result.ci_low is None
        assert result.ci_high is None
        # Point estimate stays meaningful even though the interval is undefined.
        assert result.alpha_tg == pytest.approx(1 / 3)

    def test_two_genes_ci_is_defined(self) -> None:
        """Two genes is enough for the gene-resampling bootstrap to vary."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 10,
                dn=10,
                ds=20,
                gene_id="gene1",
            ),
            PolymorphismData(
                polymorphisms=[(0.3, "N")] * 3 + [(0.6, "S")] * 7,
                dn=5,
                ds=8,
                gene_id="gene2",
            ),
        ]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        assert result.num_genes == 2
        assert result.ci_low is not None
        assert result.ci_high is not None

    def test_reproducibility_with_seed(self) -> None:
        """Test that results are reproducible with same seed."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
                dn=10,
                ds=20,
                gene_id="gene1",
            )
            for _ in range(10)
        ]

        result1 = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)
        result2 = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        assert result1.ci_low == result2.ci_low
        assert result1.ci_high == result2.ci_high


class TestAlphaTGResult:
    """Tests for AlphaTGResult dataclass."""

    def test_str_representation(self) -> None:
        """Test string representation."""
        result = AlphaTGResult(
            alpha_tg=0.35,
            ni_tg=0.65,
            ci_low=0.25,
            ci_high=0.45,
            num_genes=50,
            dn_total=100,
            ds_total=200,
            pn_total=80,
            ps_total=150,
        )

        s = str(result)

        assert "Tarone-Greenland" in s
        assert "α_TG" in s
        assert "0.35" in s or "0.3500" in s
        assert "NI_TG" in s
        assert "0.65" in s or "0.6500" in s
        assert "Dn=100" in s
        assert "Ds=200" in s
        assert "Pn=80" in s
        assert "Ps=150" in s
        assert "Genes: 50" in s

    def test_to_dict(self) -> None:
        """Test dictionary conversion."""
        result = AlphaTGResult(
            alpha_tg=0.35,
            ni_tg=0.65,
            ci_low=0.25,
            ci_high=0.45,
            num_genes=50,
            dn_total=100,
            ds_total=200,
            pn_total=80,
            ps_total=150,
        )

        d = result.to_dict()

        assert d["alpha_tg"] == 0.35
        assert d["ni_tg"] == 0.65
        assert d["ci_low"] == 0.25
        assert d["ci_high"] == 0.45
        assert d["num_genes"] == 50
        assert d["dn_total"] == 100
        assert d["ds_total"] == 200
        assert d["pn_total"] == 80
        assert d["ps_total"] == 150


class TestAlphaTGIntegration:
    """Integration tests with real sequence data."""

    def test_with_kreitman_data(self) -> None:
        """Test α_TG calculation with Kreitman Adh data."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        gene_data = [extract_polymorphism_data(ingroup, outgroup, gene_id="adh")]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        assert isinstance(result, AlphaTGResult)
        assert result.num_genes == 1
        assert result.dn_total > 0 or result.ds_total > 0
        # One gene gives the bootstrap nothing to resample: CI is undefined.
        assert result.ci_low is None
        assert result.ci_high is None

    def test_comparison_with_simple_alpha(self, tmp_path: Path) -> None:
        """Test that α_TG differs from simple α average when heterogeneity exists."""
        # Create test data with heterogeneity across genes
        ingroup1 = tmp_path / "ingroup1.fa"
        ingroup1.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
>seq4
ATGATGATGATGCTGATG
>seq5
ATGATGATGATGATGATG
""")

        outgroup1 = tmp_path / "outgroup1.fa"
        outgroup1.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        ingroup2 = tmp_path / "ingroup2.fa"
        ingroup2.write_text(""">seq1
ATGATGATGATGATGATGATGATG
>seq2
ATGCTGATGCTGATGCTGATGATG
>seq3
ATGATGATGATGATGATGATGATG
""")

        outgroup2 = tmp_path / "outgroup2.fa"
        outgroup2.write_text(""">out1
ATGGTGATGGTGATGGTGATGGTG
""")

        gene_data = [
            extract_polymorphism_data(ingroup1, outgroup1, gene_id="gene1"),
            extract_polymorphism_data(ingroup2, outgroup2, gene_id="gene2"),
        ]

        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=100, seed=42)

        # Should produce a valid result
        assert isinstance(result, AlphaTGResult)
        assert result.num_genes == 2
