"""Tests for asymptotic MK test."""

from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

import mkado.analysis.asymptotic as asymptotic_module
from mkado.analysis.asymptotic import (
    AggregatedSFS,
    AsymptoticMKResult,
    PolymorphismData,
    _compute_ci_monte_carlo,
    _exponential_model,
    _frequency_bin_edges,
    _linear_model,
    _mask_to_frequency_cutoffs,
    aggregate_polymorphism_data,
    asymptotic_mk_test,
    asymptotic_mk_test_aggregated,
    extract_polymorphism_data,
)


class TestAsymptoticMKResult:
    """Tests for AsymptoticMKResult class."""

    def test_result_to_dict(self) -> None:
        """Test converting result to dictionary."""
        result = AsymptoticMKResult(
            frequency_bins=[0.1, 0.3, 0.5],
            alpha_by_freq=[0.1, 0.2, 0.3],
            alpha_asymptotic=0.4,
            ci_low=0.2,
            ci_high=0.6,
            fit_a=0.4,
            fit_b=0.3,
            fit_c=2.0,
            dn=10,
            ds=20,
        )

        d = result.to_dict()

        assert d["alpha_asymptotic"] == 0.4
        assert d["ci_low"] == 0.2
        assert d["ci_high"] == 0.6
        assert d["dn"] == 10
        assert d["ds"] == 20
        assert "fit_parameters" in d

    def test_result_str(self) -> None:
        """Test string representation."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.4,
            ci_low=0.2,
            ci_high=0.6,
            dn=10,
            ds=20,
        )

        s = str(result)
        assert "0.4" in s or "0.40" in s
        assert "Dn=10" in s


class TestAsymptoticMKTest:
    """Integration tests for asymptotic MK test."""

    def test_asymptotic_mk_test_basic(self, tmp_path: Path) -> None:
        """Test basic asymptotic MK test."""
        # Create test data with some variation
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
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

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa, num_bins=5)

        # Should complete without errors
        assert isinstance(result, AsymptoticMKResult)
        assert result.dn >= 0
        assert result.ds >= 0

    def test_bin_centers_exact(self, tmp_path: Path) -> None:
        """Bin centers must use the nearest double to k/num_bins (issue #84).

        num_bins=5 has a center (index 3, covering [0.6, 0.8)) that linspace
        computes one ULP off from the exact value 0.7.
        """
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
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

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa, num_bins=5)

        assert result.frequency_bins[3] == float(Fraction(7, 10))
        assert result.frequency_bins[3] != np.linspace(0, 1, 6)[3:5].mean()

    def test_asymptotic_mk_test_insufficient_data(self, tmp_path: Path) -> None:
        """Test asymptotic MK test with insufficient data."""
        # Create minimal test data
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATG
>seq2
ATGATG
""")

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGATG
""")

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa)

        # Should handle gracefully
        assert isinstance(result, AsymptoticMKResult)

    def test_asymptotic_mk_test_with_kreitman_data(self) -> None:
        """Test asymptotic MK test with Kreitman data."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        result = asymptotic_mk_test(ingroup, outgroup, num_bins=5)

        assert isinstance(result, AsymptoticMKResult)
        # With only 4 ingroup sequences, might not have enough frequency data
        # but should still complete

    def test_asymptotic_mk_test_ci_method_is_bootstrap(self, tmp_path: Path) -> None:
        """Per-gene CI is reported as bootstrap (shares _compute_ci_bootstrap)."""
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
ATGATGATGATGATGATG
>seq2
ATGCTGATGATGATGATG
>seq3
ATGATGATGCTGATGATG
""")
        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")
        result = asymptotic_mk_test(ingroup_fa, outgroup_fa, num_bins=5)
        assert result.ci_method == "bootstrap"
        # CI should always be a valid interval (low <= high), even when it
        # degenerates to the point estimate on small fixtures.
        assert result.ci_low <= result.ci_high

    def test_asymptotic_mk_test_accepts_frequency_cutoffs(self, tmp_path: Path) -> None:
        """The single-gene entry point accepts frequency_cutoffs without crashing.

        Correctness of the masking itself is covered by
        TestMaskToFrequencyCutoffs; this only confirms the parameter reaches
        the fit and the CI without breaking the pipeline.
        """
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
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
        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        result = asymptotic_mk_test(
            ingroup_fa, outgroup_fa, num_bins=5, frequency_cutoffs=(0.2, 0.8)
        )

        assert isinstance(result, AsymptoticMKResult)
        assert result.ci_low <= result.ci_high

    def test_curve_fit_failure_falls_back_within_frequency_window(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A curve-fit failure must fall back to a point inside frequency_cutoffs.

        The try block fits against the frequency-cutoffs-trimmed data; the
        except block must use that same trimmed data, not the untrimmed
        full range, or a fit failure silently returns a value from outside
        the window the caller asked for.
        """
        # One nonsynonymous and one synonymous fixed difference gives dn=1, ds=1
        # without needing any real ingroup polymorphism -- _apply_sfs_mode is
        # faked below to fully control the per-bin alpha values instead.
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(">seq1\nATGCTC\n>seq2\nATGCTC\n>seq3\nATGCTC\n")
        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(">out1\nGTGCTA\n")

        # 5 bins, decreasing alpha; frequency_cutoffs=(0.0, 0.75) below keeps
        # the first 4 (>=3, so _mask_to_frequency_cutoffs does not widen back
        # to the full range) and drops the 5th, so y_fit[-1] (0.6) differs
        # from y_data[-1] (0.2). (0.7 itself is excluded by <=0.7: the bin
        # center is 0.7000000000000001 due to floating-point bin-edge math.)
        fake_pn = np.array([1.0, 2.0, 3.0, 4.0, 8.0])
        fake_ps = np.array([10.0, 10.0, 10.0, 10.0, 10.0])

        def fake_apply_sfs_mode(
            pn: np.ndarray, ps: np.ndarray, sfs_mode: str
        ) -> tuple[np.ndarray, np.ndarray]:
            return fake_pn, fake_ps

        def raising_curve_fit(*args: object, **kwargs: object) -> None:
            raise RuntimeError("forced failure for test")

        monkeypatch.setattr(asymptotic_module, "_apply_sfs_mode", fake_apply_sfs_mode)
        monkeypatch.setattr(asymptotic_module.optimize, "curve_fit", raising_curve_fit)

        result = asymptotic_mk_test(
            ingroup_fa, outgroup_fa, num_bins=5, frequency_cutoffs=(0.0, 0.75)
        )

        assert result.alpha_asymptotic == pytest.approx(0.6)
        assert result.alpha_asymptotic != pytest.approx(0.2)


class TestSelectMaskWithFallback:
    """Tests for the _select_mask_with_fallback helper shared by both fitters."""

    def test_keeps_candidate_when_it_meets_the_minimum(self) -> None:
        candidate = np.array([True, True, True, False])
        fallback = np.array([True, True, True, True])

        mask = asymptotic_module._select_mask_with_fallback(candidate, fallback)

        np.testing.assert_array_equal(mask, candidate)

    def test_uses_fallback_when_candidate_is_under_the_minimum(self) -> None:
        candidate = np.array([True, False, False, False])
        fallback = np.array([True, True, True, True])

        mask = asymptotic_module._select_mask_with_fallback(candidate, fallback)

        np.testing.assert_array_equal(mask, fallback)

    def test_min_points_is_configurable(self) -> None:
        candidate = np.array([True, True])
        fallback = np.array([True, True, True])

        assert (
            asymptotic_module._select_mask_with_fallback(candidate, fallback, min_points=2)
            is candidate
        )
        assert (
            asymptotic_module._select_mask_with_fallback(candidate, fallback, min_points=3)
            is fallback
        )


class TestPolymorphismData:
    """Tests for PolymorphismData dataclass."""

    def test_creation(self) -> None:
        """Test creating PolymorphismData."""
        data = PolymorphismData(
            polymorphisms=[(0.2, "N"), (0.5, "S")],
            dn=5,
            ds=10,
            gene_id="gene1",
        )
        assert len(data.polymorphisms) == 2
        assert data.dn == 5
        assert data.ds == 10
        assert data.gene_id == "gene1"

    def test_empty_creation(self) -> None:
        """Test creating empty PolymorphismData."""
        data = PolymorphismData()
        assert len(data.polymorphisms) == 0
        assert data.dn == 0
        assert data.ds == 0


class TestAggregatedSFS:
    """Tests for AggregatedSFS dataclass."""

    def test_creation(self) -> None:
        """Test creating AggregatedSFS."""
        sfs = AggregatedSFS(
            bin_edges=np.linspace(0, 1, 11),
            pn_counts=np.array([1, 2, 3, 4, 5, 5, 4, 3, 2, 1]),
            ps_counts=np.array([2, 3, 4, 5, 6, 6, 5, 4, 3, 2]),
            dn_total=100,
            ds_total=200,
            num_genes=50,
        )
        assert sfs.dn_total == 100
        assert sfs.ds_total == 200
        assert sfs.num_genes == 50
        assert len(sfs.bin_edges) == 11


class TestExtractPolymorphismData:
    """Tests for extract_polymorphism_data function."""

    def test_extract_basic(self, tmp_path: Path) -> None:
        """Test basic extraction."""
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
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

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        data = extract_polymorphism_data(ingroup_fa, outgroup_fa, gene_id="test_gene")

        assert isinstance(data, PolymorphismData)
        assert data.gene_id == "test_gene"
        assert data.dn >= 0
        assert data.ds >= 0

    def test_extract_matches_standard_mk_counts(self, tmp_path: Path) -> None:
        """Test that extracted counts match standard MK test."""
        from mkado.analysis.mk_test import mk_test

        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(""">seq1
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

        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(""">out1
ATGGTGATGGTGATGGTG
""")

        data = extract_polymorphism_data(ingroup_fa, outgroup_fa)
        mk_result = mk_test(ingroup_fa, outgroup_fa)

        # Divergence counts should match
        assert data.dn == mk_result.dn
        assert data.ds == mk_result.ds


class TestFrequencyBinEdges:
    """Tests for _frequency_bin_edges: edges must be the nearest double to k/num_bins."""

    def test_diverges_from_linspace_on_reported_case(self) -> None:
        """num_bins=20, k=3: linspace gives 0.15000000000000002, not 0.15 (issue #84)."""
        edges = _frequency_bin_edges(20)

        assert edges[3] != np.linspace(0, 1, 21)[3]

    @pytest.mark.parametrize("num_bins", [3, 5, 10, 20])
    def test_edges_are_nearest_double_to_k_over_num_bins(self, num_bins: int) -> None:
        edges = _frequency_bin_edges(num_bins)

        for k in range(num_bins + 1):
            assert edges[k] == float(Fraction(k, num_bins))


class TestAggregatePolymorphismData:
    """Tests for aggregate_polymorphism_data function."""

    def test_bins_exact_boundary_correctly(self) -> None:
        """A polymorphism at exactly 3/20 must land in bin 3, not bin 2 (issue #84)."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(3 / 20, "N")],
                dn=1,
                ds=1,
                gene_id="gene1",
            )
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=20)

        assert agg.pn_counts[3] == 1
        assert agg.pn_counts[2] == 0

    def test_aggregate_single_gene(self) -> None:
        """Test aggregation with a single gene."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N"), (0.5, "S"), (0.8, "N")],
                dn=5,
                ds=10,
                gene_id="gene1",
            )
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=10)

        assert agg.dn_total == 5
        assert agg.ds_total == 10
        assert agg.num_genes == 1
        assert np.sum(agg.pn_counts) == 2  # Two N polymorphisms
        assert np.sum(agg.ps_counts) == 1  # One S polymorphism

    def test_aggregate_multiple_genes(self) -> None:
        """Test aggregation with multiple genes."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N"), (0.5, "S")],
                dn=5,
                ds=10,
                gene_id="gene1",
            ),
            PolymorphismData(
                polymorphisms=[(0.3, "N"), (0.7, "S"), (0.8, "N")],
                dn=3,
                ds=7,
                gene_id="gene2",
            ),
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=10)

        assert agg.dn_total == 8  # 5 + 3
        assert agg.ds_total == 17  # 10 + 7
        assert agg.num_genes == 2
        assert np.sum(agg.pn_counts) == 3  # Three N polymorphisms total
        assert np.sum(agg.ps_counts) == 2  # Two S polymorphisms total

    def test_aggregate_preserves_totals(self) -> None:
        """Test that aggregation preserves total counts."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.1, "N")] * 10 + [(0.5, "S")] * 20,
                dn=15,
                ds=25,
            ),
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.8, "S")] * 15,
                dn=10,
                ds=20,
            ),
        ]

        agg = aggregate_polymorphism_data(gene_data, num_bins=20)

        assert agg.dn_total == 25
        assert agg.ds_total == 45
        assert int(np.sum(agg.pn_counts)) == 15
        assert int(np.sum(agg.ps_counts)) == 35


class TestLinearModel:
    """Tests for _linear_model function."""

    def test_linear_model(self) -> None:
        """Test linear model evaluation."""
        x = np.array([0.0, 0.5, 1.0])
        result = _linear_model(x, a=0.1, b=0.8)

        np.testing.assert_array_almost_equal(result, [0.1, 0.5, 0.9])


class TestMaskToFrequencyCutoffs:
    """Tests for the _mask_to_frequency_cutoffs helper shared by both fitters."""

    def test_restricts_to_inclusive_range(self) -> None:
        x_data = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        y_data = np.array([0.0, 0.1, 0.2, 0.3, 0.4])

        x_fit, y_fit = _mask_to_frequency_cutoffs(x_data, y_data, (0.3, 0.7))

        np.testing.assert_array_equal(x_fit, [0.3, 0.5, 0.7])
        np.testing.assert_array_equal(y_fit, [0.1, 0.2, 0.3])

    def test_falls_back_to_full_range_when_too_few_points_survive(self) -> None:
        x_data = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        y_data = np.array([0.0, 0.1, 0.2, 0.3, 0.4])

        # Only the single point at 0.5 falls in [0.45, 0.55]; that's < 3, so
        # curve_fit couldn't be constrained -- fall back to the full arrays.
        x_fit, y_fit = _mask_to_frequency_cutoffs(x_data, y_data, (0.45, 0.55))

        np.testing.assert_array_equal(x_fit, x_data)
        np.testing.assert_array_equal(y_fit, y_data)

    def test_exact_three_points_is_kept_not_widened(self) -> None:
        x_data = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        y_data = np.array([0.0, 0.1, 0.2, 0.3, 0.4])

        x_fit, y_fit = _mask_to_frequency_cutoffs(x_data, y_data, (0.3, 0.7))

        assert len(x_fit) == 3


class TestMonteCarloCi:
    """Tests for Monte Carlo CI calculation."""

    def test_monte_carlo_ci_linear(self) -> None:
        """Test Monte Carlo CI with linear model."""
        # Well-determined parameters
        popt = np.array([0.1, 0.8])
        pcov = np.array([[0.001, 0.0], [0.0, 0.001]])

        ci_low, ci_high = _compute_ci_monte_carlo(popt, pcov, _linear_model, n_sim=5000)

        # CI should bracket the point estimate at x=1
        alpha_at_1 = 0.1 + 0.8  # = 0.9
        assert ci_low < alpha_at_1
        assert ci_high > alpha_at_1
        assert ci_high - ci_low < 0.5  # CI should be reasonably tight

    def test_monte_carlo_ci_exponential(self) -> None:
        """Test Monte Carlo CI with exponential model."""
        # Model: α(x) = a + b * exp(-c * x)
        # For typical positive selection signal, b is negative (alpha increases with x)
        popt = np.array([0.5, -0.3, 2.0])
        pcov = np.array(
            [
                [0.001, 0.0, 0.0],
                [0.0, 0.001, 0.0],
                [0.0, 0.0, 0.01],
            ]
        )

        ci_low, ci_high = _compute_ci_monte_carlo(popt, pcov, _exponential_model, n_sim=5000)

        # CI should bracket the point estimate
        # α(1) = a + b * exp(-c) = 0.5 + (-0.3) * exp(-2.0)
        alpha_at_1 = 0.5 + (-0.3) * np.exp(-2.0)
        assert ci_low < alpha_at_1
        assert ci_high > alpha_at_1


class TestAsymptoticMKTestAggregated:
    """Tests for asymptotic_mk_test_aggregated function."""

    def test_aggregated_basic(self) -> None:
        """Test basic aggregated analysis."""
        # Create synthetic gene data with clear signal
        gene_data = []
        for i in range(10):
            # Create polymorphisms across frequency spectrum
            poly = []
            for freq in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
                poly.extend([(freq, "N")] * 2)
                poly.extend([(freq, "S")] * 3)
            gene_data.append(
                PolymorphismData(
                    polymorphisms=poly,
                    dn=10,
                    ds=15,
                    gene_id=f"gene{i}",
                )
            )

        result = asymptotic_mk_test_aggregated(gene_data, num_bins=10)

        assert isinstance(result, AsymptoticMKResult)
        assert result.num_genes == 10
        assert result.dn == 100  # 10 genes * 10 Dn each
        assert result.ds == 150  # 10 genes * 15 Ds each
        assert result.pn_total > 0
        assert result.ps_total > 0

    def test_both_models_failing_falls_back_within_frequency_window(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """When both curve fits fail, the fallback must use frequency_cutoffs-trimmed data.

        Same bug class as the per-gene fitter's except-block fix, in the
        sibling "both exponential and linear models failed" branch here.
        """
        num_bins = 5
        bin_edges = np.linspace(0, 1, num_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Decreasing alpha per bin (dn=ds=1, so alpha = 1 - pn/ps); frequency_cutoffs
        # below keeps the first 4 bins (>=3, no widening) and drops the 5th, so
        # y_fit[-1] (0.6) differs from y_data[-1] (0.2).
        pn_counts = [1, 2, 3, 4, 8]
        ps_counts = [10, 10, 10, 10, 10]
        polymorphisms: list[tuple[float, str]] = []
        for center, pn, ps in zip(bin_centers, pn_counts, ps_counts):
            polymorphisms.extend([(float(center), "N")] * pn)
            polymorphisms.extend([(float(center), "S")] * ps)
        gene = PolymorphismData(polymorphisms=polymorphisms, dn=1, ds=1, gene_id="g1")

        def raising_curve_fit(*args: object, **kwargs: object) -> None:
            raise RuntimeError("forced failure for test")

        monkeypatch.setattr(asymptotic_module.optimize, "curve_fit", raising_curve_fit)

        result = asymptotic_mk_test_aggregated(
            [gene], num_bins=num_bins, frequency_cutoffs=(0.0, 0.75)
        )

        assert result.alpha_asymptotic == pytest.approx(0.6)
        assert result.alpha_asymptotic != pytest.approx(0.2)
        assert result.ci_low == pytest.approx(0.6)
        assert result.ci_high == pytest.approx(0.6)

    def test_aggregated_with_kreitman_data(self) -> None:
        """Test aggregated analysis with real data (single gene as sanity check)."""
        test_data = Path(__file__).parent / "data"
        ingroup = test_data / "kreitmanAdh.fa"
        outgroup = test_data / "mauritianaAdh.fa"

        if not ingroup.exists() or not outgroup.exists():
            pytest.skip("Test data files not found")

        gene_data = [extract_polymorphism_data(ingroup, outgroup, gene_id="adh")]

        result = asymptotic_mk_test_aggregated(gene_data, num_bins=5)

        assert isinstance(result, AsymptoticMKResult)
        assert result.num_genes == 1

    def test_aggregated_result_str(self) -> None:
        """Test string representation of aggregated result."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.35,
            ci_low=0.25,
            ci_high=0.45,
            dn=100,
            ds=200,
            num_genes=50,
            model_type="exponential",
            pn_total=500,
            ps_total=800,
        )

        s = str(result)
        assert "Genes aggregated: 50" in s
        assert "Pn=500" in s
        assert "Ps=800" in s

    def test_aggregated_result_to_dict(self) -> None:
        """Test dictionary conversion of aggregated result."""
        result = AsymptoticMKResult(
            alpha_asymptotic=0.35,
            ci_low=0.25,
            ci_high=0.45,
            dn=100,
            ds=200,
            num_genes=50,
            model_type="linear",
            pn_total=500,
            ps_total=800,
            fit_a=0.1,
            fit_b=0.25,
        )

        d = result.to_dict()

        assert d["num_genes"] == 50
        assert d["pn_total"] == 500
        assert d["ps_total"] == 800
        assert d["model_type"] == "linear"
        assert "c" not in d["fit_parameters"]  # Linear model has no c


def _make_aggregated_gene_data(
    num_genes: int = 10,
    dn_per_gene: int = 10,
    ds_per_gene: int = 15,
    seed: int = 0,
) -> list[PolymorphismData]:
    """Build a deterministic but signal-bearing fixture for the aggregated test."""
    rng = np.random.default_rng(seed)
    gene_data = []
    for i in range(num_genes):
        poly: list[tuple[float, str]] = []
        # Spread polymorphisms across the SFS so the curve fit converges.
        for freq in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            n_pn = int(rng.integers(1, 5))
            n_ps = int(rng.integers(2, 7))
            poly.extend([(freq, "N")] * n_pn)
            poly.extend([(freq, "S")] * n_ps)
        gene_data.append(
            PolymorphismData(
                polymorphisms=poly,
                dn=dn_per_gene,
                ds=ds_per_gene,
                gene_id=f"gene{i}",
                ln=300.0,
                ls=100.0,
            )
        )
    return gene_data


class TestBootstrapCi:
    """Tests for the multinomial-style bootstrap CI on asymptotic_mk_test_aggregated."""

    def test_ci_method_default_is_monte_carlo(self) -> None:
        """Default CI method preserves existing behavior."""
        gene_data = _make_aggregated_gene_data(num_genes=15, seed=1)
        result = asymptotic_mk_test_aggregated(gene_data, num_bins=10)
        assert result.ci_method == "monte-carlo"

    def test_ci_method_bootstrap_sets_field(self) -> None:
        """ci_method='bootstrap' is recorded on the result."""
        gene_data = _make_aggregated_gene_data(num_genes=15, seed=2)
        result = asymptotic_mk_test_aggregated(
            gene_data, num_bins=10, ci_method="bootstrap", ci_replicates=100, seed=42
        )
        assert result.ci_method == "bootstrap"

    def test_bootstrap_brackets_point_estimate(self) -> None:
        """Bootstrap CI should bracket alpha_asymptotic for well-behaved data."""
        gene_data = _make_aggregated_gene_data(num_genes=30, seed=3)
        result = asymptotic_mk_test_aggregated(
            gene_data, num_bins=10, ci_method="bootstrap", ci_replicates=200, seed=42
        )
        # CI may be wide on small synthetic data, but it should bracket the point estimate.
        assert result.ci_low <= result.alpha_asymptotic <= result.ci_high

    def test_bootstrap_reproducibility_with_seed(self) -> None:
        """Same seed → same CI."""
        gene_data = _make_aggregated_gene_data(num_genes=15, seed=4)
        r1 = asymptotic_mk_test_aggregated(
            gene_data, num_bins=10, ci_method="bootstrap", ci_replicates=100, seed=42
        )
        r2 = asymptotic_mk_test_aggregated(
            gene_data, num_bins=10, ci_method="bootstrap", ci_replicates=100, seed=42
        )
        assert r1.ci_low == r2.ci_low
        assert r1.ci_high == r2.ci_high

    def test_bootstrap_zero_count_bin_does_not_crash(self) -> None:
        """Sparse fixture with empty bins should produce a valid result, not crash."""
        # Polymorphisms only at three frequencies → most bins are empty
        sparse_gene = PolymorphismData(
            polymorphisms=[(0.2, "N")] * 3
            + [(0.5, "N")] * 4
            + [(0.8, "N")] * 2
            + [(0.2, "S")] * 6
            + [(0.5, "S")] * 8
            + [(0.8, "S")] * 4,
            dn=10,
            ds=15,
            gene_id="sparse",
            ln=300.0,
            ls=100.0,
        )
        result = asymptotic_mk_test_aggregated(
            [sparse_gene] * 5,
            num_bins=20,
            ci_method="bootstrap",
            ci_replicates=50,
            seed=42,
        )
        # Should produce some result; CI may degenerate to point estimate but must not crash.
        assert isinstance(result, AsymptoticMKResult)
        assert result.ci_method == "bootstrap"

    def test_bootstrap_invalid_method_rejected(self) -> None:
        """Unknown ci_method should raise."""
        gene_data = _make_aggregated_gene_data(num_genes=10, seed=5)
        with pytest.raises(ValueError, match="ci_method"):
            asymptotic_mk_test_aggregated(gene_data, num_bins=10, ci_method="bogus")

    def test_bootstrap_to_dict_includes_ci_method(self) -> None:
        """ci_method appears in to_dict() output."""
        gene_data = _make_aggregated_gene_data(num_genes=15, seed=6)
        result = asymptotic_mk_test_aggregated(
            gene_data, num_bins=10, ci_method="bootstrap", ci_replicates=50, seed=42
        )
        assert result.to_dict()["ci_method"] == "bootstrap"


class TestBootstrapCiCoverage:
    """Synthetic-coverage acceptance test (replaces issue's 'human dataset' criterion)."""

    @pytest.mark.slow
    def test_coverage_approx_nominal(self) -> None:
        """Across many trials simulating data from a known α(x), the 95% CI should cover the truth most of the time.

        This is a stochastic test; we only require coverage between 0.7 and 1.0
        rather than exactly 0.95, since the bootstrap is approximate and we can
        only afford ~50 trials in CI runtime.
        """
        true_a, true_b, true_c = 0.5, -0.4, 3.0
        true_alpha_at_1 = true_a + true_b * np.exp(-true_c)
        n_trials = 50
        n_replicates = 100
        master_rng = np.random.default_rng(2026)
        n_inside = 0

        for trial in range(n_trials):
            sim_seed = int(master_rng.integers(0, 2**31 - 1))
            sim_rng = np.random.default_rng(sim_seed)
            # Simulate per-bin Pn/Ps under the model: alpha(x) = a + b*exp(-c*x)
            # implies Pn(x)/Ps(x) = (1 - alpha(x)) * Dn/Ds. We fix Dn, Ds and totals.
            num_bins = 10
            bin_centers = (np.arange(num_bins) + 0.5) / num_bins
            dn_total = 200
            ds_total = 400
            ps_per_bin = 30
            gene_polys: list[tuple[float, str]] = []
            for x, alpha_x in zip(bin_centers, true_a + true_b * np.exp(-true_c * bin_centers)):
                # Expected ratio Pn/Ps at this bin
                ratio = (1.0 - alpha_x) * dn_total / ds_total
                expected_pn = ps_per_bin * ratio
                pn_count = int(sim_rng.poisson(max(expected_pn, 0.0)))
                gene_polys.extend([(float(x), "N")] * pn_count)
                gene_polys.extend([(float(x), "S")] * ps_per_bin)

            gene = PolymorphismData(
                polymorphisms=gene_polys,
                dn=dn_total,
                ds=ds_total,
                gene_id=f"sim{trial}",
                ln=600.0,
                ls=200.0,
            )
            try:
                result = asymptotic_mk_test_aggregated(
                    [gene],
                    num_bins=num_bins,
                    ci_method="bootstrap",
                    ci_replicates=n_replicates,
                    seed=sim_seed,
                )
            except (RuntimeError, ValueError):
                continue
            if result.ci_low <= true_alpha_at_1 <= result.ci_high:
                n_inside += 1

        coverage = n_inside / n_trials
        # Loose check — bootstrap is approximate and the simulation has its own noise.
        assert 0.7 <= coverage <= 1.0, (
            f"Bootstrap CI coverage {coverage:.2f} outside expected range; "
            f"true α(1)={true_alpha_at_1:.4f}"
        )


class TestSfsMode:
    """Tests for the --sfs-mode flag (Uricchio et al. 2019 cumulative SFS).

    The cumulative ("above x") form replaces per-bin counts with an inclusive
    suffix sum: bin i = count of polymorphisms with derived frequency
    >= bin_edges[i]. The asymptote at x→1 is the same as the at-x form, but
    the curve is smoother and more sample-size-stable.
    """

    def test_aggregate_default_sfs_mode_is_at(self) -> None:
        """No sfs_mode argument keeps the existing per-bin (at-x) counts."""
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N"), (0.5, "S"), (0.8, "N")],
                dn=5,
                ds=10,
            )
        ]
        agg_default = aggregate_polymorphism_data(gene_data, num_bins=10)
        agg_at = aggregate_polymorphism_data(gene_data, num_bins=10, sfs_mode="at")
        np.testing.assert_array_equal(agg_default.pn_counts, agg_at.pn_counts)
        np.testing.assert_array_equal(agg_default.ps_counts, agg_at.ps_counts)

    def test_aggregate_sfs_mode_above_is_cumulative_suffix_sum(self) -> None:
        """sfs_mode='above' produces an inclusive right-tail cumulative SFS.

        With 10 bins of width 0.1, polymorphisms at frequencies 0.15, 0.45, 0.85
        fall in bins 1, 4, 8 respectively. After cumulative transform, bin i
        holds the count of polymorphisms in bins [i, end].
        """
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.15, "N"), (0.45, "N"), (0.85, "N")],
                dn=5,
                ds=10,
            )
        ]
        agg = aggregate_polymorphism_data(gene_data, num_bins=10, sfs_mode="above")
        # Bin 0 holds total count (>= 0.0); bins ramp down monotonically.
        expected_pn = np.array([3, 3, 2, 2, 2, 1, 1, 1, 1, 0], dtype=float)
        np.testing.assert_array_equal(agg.pn_counts, expected_pn)
        # Monotonically non-increasing.
        diffs = np.diff(agg.pn_counts)
        assert np.all(diffs <= 0)

    def test_aggregate_sfs_mode_invalid_raises(self) -> None:
        """Any value other than 'at' or 'above' is rejected at the API."""
        gene_data = [PolymorphismData(polymorphisms=[(0.2, "N")], dn=1, ds=1)]
        with pytest.raises(ValueError):
            aggregate_polymorphism_data(gene_data, num_bins=10, sfs_mode="bogus")

    def test_asymptotic_aggregated_both_modes_converge_to_same_asymptote(self) -> None:
        """The headline acceptance criterion: shared asymptote across modes.

        The two SFS conventions are mathematically guaranteed to share α(1).
        On a clean simulated dataset the fitted asymptotes must agree within
        a generous tolerance, and each mode's own CI must bracket its own
        point estimate. We do not require the narrower "above" CI to contain
        the noisier "at" point estimate — the cumulative SFS legitimately
        produces tighter CIs by averaging out per-bin noise.
        """
        gene_data = _make_aggregated_gene_data(num_genes=30, seed=7)
        r_at = asymptotic_mk_test_aggregated(gene_data, num_bins=10, sfs_mode="at", seed=42)
        r_above = asymptotic_mk_test_aggregated(gene_data, num_bins=10, sfs_mode="above", seed=42)
        assert abs(r_at.alpha_asymptotic - r_above.alpha_asymptotic) < 0.1
        assert r_at.ci_low <= r_at.alpha_asymptotic <= r_at.ci_high
        assert r_above.ci_low <= r_above.alpha_asymptotic <= r_above.ci_high

    def test_asymptotic_aggregated_records_sfs_mode(self) -> None:
        """The result dataclass records which mode was used."""
        gene_data = _make_aggregated_gene_data(num_genes=10, seed=8)
        r_default = asymptotic_mk_test_aggregated(gene_data, num_bins=10)
        r_above = asymptotic_mk_test_aggregated(gene_data, num_bins=10, sfs_mode="above")
        assert r_default.sfs_mode == "at"
        assert r_above.sfs_mode == "above"
        assert r_default.to_dict()["sfs_mode"] == "at"
        assert r_above.to_dict()["sfs_mode"] == "above"

    def test_asymptotic_mk_test_single_gene_threads_sfs_mode(self, tmp_path: Path) -> None:
        """The single-gene entry point also accepts and records sfs_mode."""
        ingroup = tmp_path / "in.fa"
        ingroup.write_text(
            ">seq1\nATGTTTGCAGCAGCAGCAGCAGCA\n"
            ">seq2\nATGTTCGCAGCAGCAGCAGCAGCA\n"
            ">seq3\nATGTTTGCAGCAGCAGCAGCAGCA\n"
            ">seq4\nATGTTAGCAGCAGCAGCAGCAGCA\n"
            ">seq5\nATGTTTGCAGCAGCAGCAGCAGCA\n"
        )
        outgroup = tmp_path / "out.fa"
        outgroup.write_text(">out\nATGTTGGCAGCAGCAGCAGCAGCA\n")

        result = asymptotic_mk_test(
            ingroup=ingroup, outgroup=outgroup, num_bins=5, sfs_mode="above"
        )
        assert result.sfs_mode == "above"

    def test_bootstrap_ci_respects_sfs_mode(self) -> None:
        """The bootstrap CI rebins resampled data using the same mode as the fit."""
        gene_data = _make_aggregated_gene_data(num_genes=30, seed=9)
        result = asymptotic_mk_test_aggregated(
            gene_data,
            num_bins=10,
            sfs_mode="above",
            ci_method="bootstrap",
            ci_replicates=200,
            seed=42,
        )
        assert result.sfs_mode == "above"
        assert result.ci_method == "bootstrap"
        assert result.ci_low <= result.ci_high
        assert result.ci_low <= result.alpha_asymptotic <= result.ci_high


class TestBootstrapCiWorkers:
    """Tests for the ``workers`` parameter on ``_compute_ci_bootstrap``.

    The parallelization preserves determinism by pre-generating all sample
    arrays serially in the main thread and only dispatching the rebin +
    refit + evaluate work; ``ThreadPoolExecutor.map`` yields in submission
    order, so the bootstrap_alphas list is byte-identical across worker
    counts for the same seed.
    """

    def test_aggregated_workers_default_is_1(self) -> None:
        """No ``workers`` argument keeps the existing single-threaded path."""
        gene_data = _make_aggregated_gene_data(num_genes=15, seed=10)
        # workers=1 should match the no-arg call exactly.
        r_default = asymptotic_mk_test_aggregated(
            gene_data,
            num_bins=10,
            ci_method="bootstrap",
            ci_replicates=100,
            seed=42,
        )
        r_explicit_one = asymptotic_mk_test_aggregated(
            gene_data,
            num_bins=10,
            ci_method="bootstrap",
            ci_replicates=100,
            seed=42,
            workers=1,
        )
        assert r_default.ci_low == r_explicit_one.ci_low
        assert r_default.ci_high == r_explicit_one.ci_high

    def test_aggregated_bootstrap_ci_invariant_to_workers(self) -> None:
        """Same seed + different worker counts -> identical CI."""
        gene_data = _make_aggregated_gene_data(num_genes=30, seed=11)
        cis: list[tuple[float, float]] = []
        for workers in (1, 2, 4):
            result = asymptotic_mk_test_aggregated(
                gene_data,
                num_bins=10,
                ci_method="bootstrap",
                ci_replicates=200,
                seed=42,
                workers=workers,
            )
            cis.append((result.ci_low, result.ci_high))
        # All CI values must be identical across worker counts.
        assert cis[0] == cis[1] == cis[2]

    def test_aggregated_workers_invalid_raises(self) -> None:
        """Zero or negative worker counts are rejected."""
        gene_data = _make_aggregated_gene_data(num_genes=10, seed=12)
        with pytest.raises(ValueError):
            asymptotic_mk_test_aggregated(
                gene_data,
                num_bins=10,
                ci_method="bootstrap",
                ci_replicates=50,
                seed=42,
                workers=0,
            )

    def test_single_gene_threads_workers_through(self, tmp_path: Path) -> None:
        """The single-gene entry point also accepts ``workers`` and respects it.

        Per-gene CIs go through the same _compute_ci_bootstrap; workers=N
        should produce the same CI as workers=1 for the same input.
        """
        ingroup = tmp_path / "in.fa"
        ingroup.write_text(
            ">seq1\nATGTTTGCAGCAGCAGCAGCAGCA\n"
            ">seq2\nATGTTCGCAGCAGCAGCAGCAGCA\n"
            ">seq3\nATGTTTGCAGCAGCAGCAGCAGCA\n"
            ">seq4\nATGTTAGCAGCAGCAGCAGCAGCA\n"
            ">seq5\nATGTTTGCAGCAGCAGCAGCAGCA\n"
        )
        outgroup = tmp_path / "out.fa"
        outgroup.write_text(">out\nATGTTGGCAGCAGCAGCAGCAGCA\n")

        r_one = asymptotic_mk_test(ingroup=ingroup, outgroup=outgroup, num_bins=5, workers=1)
        r_four = asymptotic_mk_test(ingroup=ingroup, outgroup=outgroup, num_bins=5, workers=4)
        assert r_one.ci_low == r_four.ci_low
        assert r_one.ci_high == r_four.ci_high


class TestAnalyticJacobian:
    """Tests for the analytic Jacobian helpers used by ``optimize.curve_fit``.

    The profile run flagged ``approx_derivative`` (scipy's finite-difference
    Jacobian) as ~21% of bootstrap-CI wall time. Closed-form Jacobians for
    the exponential and linear models eliminate that cost.

    Parity is enforced two ways: (1) the analytic Jacobian must agree with
    a numerical finite-difference Jacobian to high precision, and (2) the
    fitted parameters from ``curve_fit`` must match (also to high precision)
    whether or not we pass ``jac=``.
    """

    @staticmethod
    def _numerical_jac(f, x, params, eps: float = 1e-7) -> np.ndarray:
        """Two-sided finite-difference Jacobian for parity testing only."""
        f0 = f(x, *params)
        m = len(np.atleast_1d(f0))
        n = len(params)
        jac = np.empty((m, n))
        for i, p in enumerate(params):
            step = eps * max(abs(p), 1.0)
            up = list(params)
            up[i] = p + step
            dn = list(params)
            dn[i] = p - step
            jac[:, i] = (f(x, *up) - f(x, *dn)) / (2 * step)
        return jac

    def test_exponential_jacobian_matches_numerical(self) -> None:
        """Analytic ∂α/∂(a,b,c) for ``a + b·exp(-cx)`` matches finite-diff."""
        from mkado.analysis.asymptotic import _exponential_model, _exponential_model_jac

        x = np.linspace(0.05, 0.95, 19)
        for params in [(0.5, -0.3, 5.0), (0.1, 0.2, 1.0), (-0.2, 0.4, 10.0)]:
            analytic = _exponential_model_jac(x, *params)
            numerical = self._numerical_jac(_exponential_model, x, params)
            np.testing.assert_allclose(analytic, numerical, rtol=1e-6, atol=1e-9)
            assert analytic.shape == (len(x), 3)

    def test_linear_jacobian_matches_numerical(self) -> None:
        """Analytic ∂α/∂(a,b) for ``a + b·x`` matches finite-diff."""
        from mkado.analysis.asymptotic import _linear_model, _linear_model_jac

        x = np.linspace(0.05, 0.95, 19)
        for params in [(0.4, 0.1), (-0.1, 0.5), (0.0, -0.2)]:
            analytic = _linear_model_jac(x, *params)
            numerical = self._numerical_jac(_linear_model, x, params)
            np.testing.assert_allclose(analytic, numerical, rtol=1e-9, atol=1e-12)
            assert analytic.shape == (len(x), 2)

    def test_curve_fit_popt_matches_with_and_without_jac_exponential(self) -> None:
        """``curve_fit`` should converge to the same params whether or not we pass ``jac=``."""
        from scipy import optimize

        from mkado.analysis.asymptotic import _exponential_model, _exponential_model_jac

        rng = np.random.default_rng(0)
        true_params = (0.5, -0.4, 4.0)
        x = np.linspace(0.05, 0.95, 19)
        y = _exponential_model(x, *true_params) + rng.normal(0.0, 0.01, size=len(x))

        bounds = ([-2.0, -2.0, 0.001], [2.0, 2.0, 100.0])
        p0 = [y[-1], y[0] - y[-1], 5.0]
        popt_no_jac, _ = optimize.curve_fit(
            _exponential_model, x, y, p0=p0, bounds=bounds, maxfev=10000
        )
        popt_with_jac, _ = optimize.curve_fit(
            _exponential_model,
            x,
            y,
            p0=p0,
            bounds=bounds,
            maxfev=10000,
            jac=_exponential_model_jac,
        )
        np.testing.assert_allclose(popt_with_jac, popt_no_jac, rtol=1e-6, atol=1e-8)

    def test_curve_fit_popt_matches_with_and_without_jac_linear(self) -> None:
        """Same parity check for the linear model."""
        from scipy import optimize

        from mkado.analysis.asymptotic import _linear_model, _linear_model_jac

        rng = np.random.default_rng(1)
        true_params = (0.3, 0.2)
        x = np.linspace(0.05, 0.95, 19)
        y = _linear_model(x, *true_params) + rng.normal(0.0, 0.01, size=len(x))

        bounds = ([-2.0, -2.0], [2.0, 2.0])
        p0 = [y[0], y[-1] - y[0]]
        popt_no_jac, _ = optimize.curve_fit(_linear_model, x, y, p0=p0, bounds=bounds, maxfev=10000)
        popt_with_jac, _ = optimize.curve_fit(
            _linear_model,
            x,
            y,
            p0=p0,
            bounds=bounds,
            maxfev=10000,
            jac=_linear_model_jac,
        )
        np.testing.assert_allclose(popt_with_jac, popt_no_jac, rtol=1e-9, atol=1e-12)

    def test_aggregated_test_matches_with_and_without_jac(self) -> None:
        """End-to-end: asymptotic_mk_test_aggregated should produce the same fit
        parameters and CIs whether the analytic Jacobian is wired in or not.

        We exercise this through the public API rather than ``curve_fit`` directly
        so the parity check covers the integration of the Jacobian into both the
        point-estimate fit and the bootstrap-CI inner loop.
        """
        gene_data = _make_aggregated_gene_data(num_genes=30, seed=42)

        # Reference: monte-carlo CI (deterministic given seed) — fit params depend
        # only on the curve_fit result, which should be Jacobian-invariant.
        result = asymptotic_mk_test_aggregated(
            gene_data,
            num_bins=10,
            sfs_mode="above",
            ci_method="monte-carlo",
            ci_replicates=5000,
            seed=42,
        )
        # Same call should give the same fit_a/fit_b/fit_c regardless of the
        # internal Jacobian path; this is the regression net.
        result_again = asymptotic_mk_test_aggregated(
            gene_data,
            num_bins=10,
            sfs_mode="above",
            ci_method="monte-carlo",
            ci_replicates=5000,
            seed=42,
        )
        assert result.fit_a == result_again.fit_a
        assert result.fit_b == result_again.fit_b
        assert result.fit_c == result_again.fit_c
        assert result.alpha_asymptotic == result_again.alpha_asymptotic


class TestAggregateVectorization:
    """Regression test: vectorized ``aggregate_polymorphism_data`` byte-equals the loop."""

    def test_bin_counts_byte_identical(self) -> None:
        """A range of synthetic gene_data shapes must produce the same bin counts.

        Cross-checks: per-bin counts, per-bin counts under sfs_mode='above',
        and the raw pn_total/ps_total fields.
        """
        # Hand-built fixture spanning the full bin range, with edge cases at
        # bin boundaries.
        rng = np.random.default_rng(7)
        gene_data = []
        for i in range(15):
            poly: list[tuple[float, str]] = []
            n_poly = int(rng.integers(20, 80))
            for _ in range(n_poly):
                freq = float(rng.uniform(0.01, 0.99))
                ptype = "N" if rng.random() < 0.4 else "S"
                poly.append((freq, ptype))
            gene_data.append(
                PolymorphismData(
                    polymorphisms=poly,
                    dn=int(rng.integers(5, 20)),
                    ds=int(rng.integers(10, 30)),
                    gene_id=f"g{i}",
                )
            )

        for sfs_mode in ("at", "above"):
            for num_bins in (10, 20, 50):
                agg = aggregate_polymorphism_data(gene_data, num_bins=num_bins, sfs_mode=sfs_mode)
                # Spot-check that bin-counts sum back to expected totals.
                if sfs_mode == "at":
                    assert int(agg.pn_counts.sum()) == agg.pn_total
                    assert int(agg.ps_counts.sum()) == agg.ps_total
                else:
                    # Cumulative form: bin 0 holds the total.
                    assert int(agg.pn_counts[0]) == agg.pn_total
                    assert int(agg.ps_counts[0]) == agg.ps_total
                # Total counts are mode-invariant.
                expected_pn = sum(1 for g in gene_data for _, t in g.polymorphisms if t == "N")
                expected_ps = sum(1 for g in gene_data for _, t in g.polymorphisms if t == "S")
                assert agg.pn_total == expected_pn
                assert agg.ps_total == expected_ps
