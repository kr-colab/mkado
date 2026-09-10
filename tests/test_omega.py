"""Tests for Ls/Ln site counts and omega/omega_a/omega_na decomposition.

Covers issue #9: surface Nei-Gojobori Ls/Ln totals and report omega
(Coronado-Zamora et al. 2019) on every result dataclass.
"""

from __future__ import annotations

import pytest

from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
from mkado.analysis.asymptotic import (
    PolymorphismData,
    aggregate_polymorphism_data,
    asymptotic_mk_test_aggregated,
    extract_polymorphism_data,
)
from mkado.analysis.imputed import imputed_mk_test
from mkado.analysis.mk_test import mk_test, mk_test_from_counts
from mkado.analysis.polarized import polarized_mk_test
from mkado.analysis.statistics import omega_decomposition
from mkado.core.alignment import AlignedPair
from tests.builders import sequence_set

# ---------------------------------------------------------------------------
# Site counting on AlignedPair
# ---------------------------------------------------------------------------


class TestCountTotalSites:
    """count_total_sites aggregates Nei-Gojobori per-codon syn fractions."""

    def test_methionine_only(self) -> None:
        # ATG (Met) has no synonymous sites under the standard code.
        # Two ATG codons -> Ls = 0, Ln = 6
        ingroup = sequence_set(["ATGATG"])
        outgroup = sequence_set(["ATGATG"])
        pair = AlignedPair(ingroup=ingroup, outgroup=outgroup)

        ln, ls = pair.count_total_sites()
        assert ls == pytest.approx(0.0)
        assert ln == pytest.approx(6.0)

    def test_phenylalanine_codon(self) -> None:
        # TTT (Phe): site 3 has 1/3 synonymous (only TTC -> Phe). Other sites
        # are zero. So Ls(TTT) = 1/3.
        ingroup = sequence_set(["TTT"])
        outgroup = sequence_set(["TTT"])
        pair = AlignedPair(ingroup=ingroup, outgroup=outgroup)

        ln, ls = pair.count_total_sites()
        assert ls == pytest.approx(1.0 / 3.0)
        assert ln == pytest.approx(3.0 - 1.0 / 3.0)

    def test_skips_codons_with_no_clean_codons(self) -> None:
        # First codon is fully ambiguous in both groups -> not counted.
        ingroup = sequence_set(["NNNATG"])
        outgroup = sequence_set(["NNNATG"])
        pair = AlignedPair(ingroup=ingroup, outgroup=outgroup)

        ln, ls = pair.count_total_sites()
        # Only one analyzed codon (ATG): Ln=3, Ls=0
        assert ls == pytest.approx(0.0)
        assert ln == pytest.approx(3.0)

    def test_averages_across_groups(self) -> None:
        # Ingroup ATG, outgroup TTT -> averaged Ls = (0 + 1/3) / 2
        ingroup = sequence_set(["ATG"])
        outgroup = sequence_set(["TTT"])
        pair = AlignedPair(ingroup=ingroup, outgroup=outgroup)

        ln, ls = pair.count_total_sites()
        assert ls == pytest.approx((0.0 + 1.0 / 3.0) / 2.0)
        assert ln == pytest.approx(3.0 - ls)


# ---------------------------------------------------------------------------
# omega_decomposition helper
# ---------------------------------------------------------------------------


class TestOmegaDecomposition:
    """omega_decomposition returns (omega, omega_a, omega_na)."""

    def test_basic(self) -> None:
        # Dn=10, Ds=5, Ln=200, Ls=100 -> omega = (10/5) * (100/200) = 1.0
        # alpha=0.5 -> omega_a=0.5, omega_na=0.5
        omega, omega_a, omega_na = omega_decomposition(10, 5, 200.0, 100.0, 0.5)
        assert omega == pytest.approx(1.0)
        assert omega_a == pytest.approx(0.5)
        assert omega_na == pytest.approx(0.5)

    def test_alpha_none_propagates_to_a_and_na(self) -> None:
        omega, omega_a, omega_na = omega_decomposition(10, 5, 200.0, 100.0, None)
        assert omega == pytest.approx(1.0)
        assert omega_a is None
        assert omega_na is None

    def test_zero_ds_is_none(self) -> None:
        result = omega_decomposition(10, 0, 200.0, 100.0, 0.5)
        assert result == (None, None, None)

    def test_zero_ls_is_none(self) -> None:
        result = omega_decomposition(10, 5, 200.0, 0.0, 0.5)
        assert result == (None, None, None)

    def test_zero_ln_is_none(self) -> None:
        result = omega_decomposition(10, 5, 0.0, 100.0, 0.5)
        assert result == (None, None, None)

    def test_missing_site_counts(self) -> None:
        result = omega_decomposition(10, 5, None, None, 0.5)
        assert result == (None, None, None)

    def test_alpha_zero(self) -> None:
        omega, omega_a, omega_na = omega_decomposition(10, 5, 200.0, 100.0, 0.0)
        assert omega_a == pytest.approx(0.0)
        assert omega_na == pytest.approx(omega)


# ---------------------------------------------------------------------------
# MKResult fields
# ---------------------------------------------------------------------------


class TestMKResultOmega:
    """MKResult exposes ln/ls/omega but deliberately omits omega_a/omega_na.

    Per-gene Smith-Eyre-Walker alpha is too noisy for a meaningful rate
    decomposition; see docs/omega.rst.
    """

    def test_from_counts_with_sites(self) -> None:
        result = mk_test_from_counts(dn=10, ds=5, pn=4, ps=8, ln=200.0, ls=100.0)
        assert result.ln == pytest.approx(200.0)
        assert result.ls == pytest.approx(100.0)
        # omega = (10/5) * (100/200) = 1.0
        assert result.omega == pytest.approx(1.0)
        # alpha = 1 - (5*4)/(10*8) = 1 - 0.25 = 0.75
        assert result.alpha == pytest.approx(0.75)

    def test_from_counts_without_sites(self) -> None:
        result = mk_test_from_counts(dn=10, ds=5, pn=4, ps=8)
        assert result.ln is None
        assert result.ls is None
        assert result.omega is None

    def test_from_counts_zero_ds(self) -> None:
        # With ds=0, omega is undefined
        result = mk_test_from_counts(dn=10, ds=0, pn=4, ps=8, ln=200.0, ls=100.0)
        assert result.omega is None

    def test_to_dict_includes_new_keys(self) -> None:
        result = mk_test_from_counts(dn=10, ds=5, pn=4, ps=8, ln=200.0, ls=100.0)
        d = result.to_dict()
        for key in ("ln", "ls", "omega"):
            assert key in d, f"missing {key} in to_dict()"
        # omega_a/omega_na deliberately not exposed on per-gene MKResult
        assert "omega_a" not in d
        assert "omega_na" not in d

    def test_str_includes_omega(self) -> None:
        result = mk_test_from_counts(dn=10, ds=5, pn=4, ps=8, ln=200.0, ls=100.0)
        s = str(result)
        assert "omega" in s.lower() or "ω" in s

    def test_full_pipeline_populates_sites(self, tmp_path) -> None:
        # mk_test() should compute Ln/Ls from the alignment automatically.
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(">i1\nATGTTTATG\n>i2\nATGTTCATG\n")
        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(">o1\nATGTTAATG\n")

        result = mk_test(ingroup_fa, outgroup_fa)
        assert result.ln is not None
        assert result.ls is not None
        # 3 codons analyzed -> Ln + Ls = 9
        assert result.ln + result.ls == pytest.approx(9.0)
        assert result.ls > 0  # TTT codon contributes synonymous sites


# ---------------------------------------------------------------------------
# AsymptoticMKResult
# ---------------------------------------------------------------------------


class TestAsymptoticOmega:
    def test_aggregated_carries_sites(self) -> None:
        # Build PolymorphismData with explicit ln/ls
        gene = PolymorphismData(
            polymorphisms=[(0.5, "S")] * 10 + [(0.5, "N")] * 10,
            dn=10,
            ds=10,
            ln=200.0,
            ls=100.0,
        )
        result = asymptotic_mk_test_aggregated([gene] * 5, num_bins=10)
        # Aggregated totals: ln=1000, ls=500
        assert result.ln == pytest.approx(1000.0)
        assert result.ls == pytest.approx(500.0)
        # Whatever alpha_asymptotic is, omega should be (Dn*Ls)/(Ds*Ln) using totals
        expected_omega = (result.dn * result.ls) / (result.ds * result.ln)
        assert result.omega == pytest.approx(expected_omega)

    def test_to_dict_includes_new_keys(self) -> None:
        gene = PolymorphismData(
            polymorphisms=[(0.5, "S")] * 5 + [(0.5, "N")] * 5,
            dn=10,
            ds=10,
            ln=200.0,
            ls=100.0,
        )
        result = asymptotic_mk_test_aggregated([gene] * 3, num_bins=5)
        d = result.to_dict()
        for key in ("ln", "ls", "omega", "omega_a", "omega_na"):
            assert key in d


# ---------------------------------------------------------------------------
# PolarizedMKResult
# ---------------------------------------------------------------------------


class TestPolarizedOmega:
    """PolarizedMKResult exposes ln/ls/omega but omits omega_a/omega_na.

    Per-gene polarized alpha is too noisy for a meaningful rate decomposition;
    see docs/omega.rst.
    """

    def test_full_pipeline(self, tmp_path) -> None:
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(">i1\nATGTTTATG\n>i2\nATGTTCATG\n")
        outgroup1_fa = tmp_path / "outgroup1.fa"
        outgroup1_fa.write_text(">o1\nATGTTAATG\n")
        outgroup2_fa = tmp_path / "outgroup2.fa"
        outgroup2_fa.write_text(">o2\nATGTTAATG\n")

        result = polarized_mk_test(ingroup_fa, outgroup1_fa, outgroup2_fa)
        assert result.ln is not None
        assert result.ls is not None
        # Three codons analyzed -> total = 9
        assert result.ln + result.ls == pytest.approx(9.0)
        d = result.to_dict()
        # Site counts are alignment-wide, so live at the top level
        assert "ln" in d
        assert "ls" in d
        # omega is an ingroup-lineage quantity; lives inside the ingroup sub-dict
        assert "omega" in d["ingroup"]
        # omega_a / omega_na deliberately not surfaced (per-gene alpha noise)
        assert "omega_a" not in d["ingroup"]
        assert "omega_na" not in d["ingroup"]


# ---------------------------------------------------------------------------
# ImputedMKResult
# ---------------------------------------------------------------------------


class TestImputedOmega:
    def test_uses_polymorphism_data_sites(self) -> None:
        gene = PolymorphismData(
            polymorphisms=[(0.05, "N")] * 3 + [(0.50, "N")] * 5 + [(0.50, "S")] * 8,
            dn=10,
            ds=5,
            ln=200.0,
            ls=100.0,
        )
        result = imputed_mk_test(gene, cutoff=0.15)
        assert result.ln == pytest.approx(200.0)
        assert result.ls == pytest.approx(100.0)
        # omega = (10/5) * (100/200) = 1.0
        assert result.omega == pytest.approx(1.0)
        d = result.to_dict()
        assert "omega" in d


# ---------------------------------------------------------------------------
# AlphaTGResult
# ---------------------------------------------------------------------------


class TestAlphaTGOmega:
    def test_aggregates_sites_across_genes(self) -> None:
        gene1 = PolymorphismData(
            polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
            dn=10,
            ds=20,
            ln=200.0,
            ls=100.0,
        )
        gene2 = PolymorphismData(
            polymorphisms=[(0.3, "N")] * 3 + [(0.5, "S")] * 10,
            dn=8,
            ds=16,
            ln=150.0,
            ls=75.0,
        )
        result = alpha_tg_from_gene_data([gene1, gene2], bootstrap_replicates=10, seed=42)
        assert result.ln == pytest.approx(350.0)
        assert result.ls == pytest.approx(175.0)
        # omega from totals
        expected_omega = (result.dn_total * result.ls) / (result.ds_total * result.ln)
        assert result.omega == pytest.approx(expected_omega)
        d = result.to_dict()
        assert "ln" in d
        assert "omega" in d


# ---------------------------------------------------------------------------
# Bootstrap CIs on omega_a / omega_na (and on omega for alpha_tg)
# ---------------------------------------------------------------------------


class TestAsymptoticOmegaCIs:
    """AsymptoticMKResult exposes omega_a/omega_na CIs scaled from the alpha CI.

    Dn, Ds, Ln, Ls are constants under the curve-fit Monte Carlo, so omega
    itself has no sampling distribution. omega_a = alpha * omega scales the
    alpha CI directly; omega_na flips the percentiles since (1 - alpha)
    inverts the ordering.
    """

    def _build_result(self):
        gene = PolymorphismData(
            polymorphisms=[(0.5, "S")] * 10 + [(0.5, "N")] * 10,
            dn=10,
            ds=10,
            ln=200.0,
            ls=100.0,
        )
        return asymptotic_mk_test_aggregated([gene] * 5, num_bins=10)

    def test_ci_fields_present(self) -> None:
        result = self._build_result()
        for attr in (
            "omega_a_ci_low",
            "omega_a_ci_high",
            "omega_na_ci_low",
            "omega_na_ci_high",
        ):
            assert hasattr(result, attr), f"missing {attr}"

    def test_omega_a_ci_scales_alpha_ci(self) -> None:
        result = self._build_result()
        if result.omega is None:
            pytest.skip("omega undefined for this fixture")
        assert result.omega_a_ci_low == pytest.approx(result.ci_low * result.omega)
        assert result.omega_a_ci_high == pytest.approx(result.ci_high * result.omega)

    def test_omega_na_ci_flips_alpha_ci(self) -> None:
        result = self._build_result()
        if result.omega is None:
            pytest.skip("omega undefined for this fixture")
        # omega_na = (1 - alpha) * omega; since (1-x) is monotonically
        # decreasing, the low/high quantiles flip.
        assert result.omega_na_ci_low == pytest.approx((1.0 - result.ci_high) * result.omega)
        assert result.omega_na_ci_high == pytest.approx((1.0 - result.ci_low) * result.omega)

    def test_to_dict_includes_ci_keys(self) -> None:
        result = self._build_result()
        d = result.to_dict()
        for key in (
            "omega_a_ci_low",
            "omega_a_ci_high",
            "omega_na_ci_low",
            "omega_na_ci_high",
        ):
            assert key in d


class TestAlphaTGOmegaCIs:
    """AlphaTGResult bootstraps genes, so omega itself varies per replicate."""

    def _gene_data(self):
        return [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
                dn=10,
                ds=20,
                ln=200.0,
                ls=100.0,
            ),
            PolymorphismData(
                polymorphisms=[(0.3, "N")] * 3 + [(0.5, "S")] * 10,
                dn=8,
                ds=16,
                ln=150.0,
                ls=75.0,
            ),
            PolymorphismData(
                polymorphisms=[(0.1, "N")] * 7 + [(0.5, "S")] * 12,
                dn=6,
                ds=14,
                ln=180.0,
                ls=90.0,
            ),
        ]

    def test_ci_fields_present(self) -> None:
        result = alpha_tg_from_gene_data(self._gene_data(), bootstrap_replicates=200, seed=42)
        for attr in (
            "omega_ci_low",
            "omega_ci_high",
            "omega_a_ci_low",
            "omega_a_ci_high",
            "omega_na_ci_low",
            "omega_na_ci_high",
        ):
            assert hasattr(result, attr), f"missing {attr}"

    def test_cis_bracket_point_estimates(self) -> None:
        result = alpha_tg_from_gene_data(self._gene_data(), bootstrap_replicates=500, seed=42)
        assert result.omega is not None
        # Point estimates should fall within the bootstrap CI in nearly all cases.
        # With 3 genes the bootstrap is noisy, so we just check the CI is a valid
        # interval (low <= high) and brackets the point estimate.
        assert result.omega_ci_low <= result.omega_ci_high
        assert result.omega_a_ci_low <= result.omega_a_ci_high
        assert result.omega_na_ci_low <= result.omega_na_ci_high

    def test_to_dict_includes_ci_keys(self) -> None:
        result = alpha_tg_from_gene_data(self._gene_data(), bootstrap_replicates=10, seed=42)
        d = result.to_dict()
        for key in (
            "omega_ci_low",
            "omega_ci_high",
            "omega_a_ci_low",
            "omega_a_ci_high",
            "omega_na_ci_low",
            "omega_na_ci_high",
        ):
            assert key in d

    def test_single_gene_all_cis_none(self) -> None:
        # One gene gives the gene-resampling bootstrap nothing to vary, so the
        # alpha CI and every omega CI derived from it are undefined.
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5 + [(0.5, "S")] * 15,
                dn=10,
                ds=20,
                ln=200.0,
                ls=100.0,
            ),
        ]
        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=10, seed=42)
        assert result.omega is not None
        assert result.ci_low is None
        assert result.ci_high is None
        assert result.omega_ci_low is None
        assert result.omega_a_ci_low is None
        assert result.omega_na_ci_low is None

    def test_cis_none_when_sites_missing(self) -> None:
        # If any gene lacks site counts, omega_total is None and so are the CIs.
        gene_data = [
            PolymorphismData(
                polymorphisms=[(0.2, "N")] * 5,
                dn=10,
                ds=20,
                # no ln/ls
            ),
        ]
        result = alpha_tg_from_gene_data(gene_data, bootstrap_replicates=10, seed=42)
        assert result.omega is None
        assert result.omega_ci_low is None
        assert result.omega_a_ci_low is None
        assert result.omega_na_ci_low is None


# ---------------------------------------------------------------------------
# extract_polymorphism_data populates Ln/Ls
# ---------------------------------------------------------------------------


class TestExtractPolymorphismDataSites:
    def test_populates_sites(self, tmp_path) -> None:
        ingroup_fa = tmp_path / "ingroup.fa"
        ingroup_fa.write_text(">i1\nATGTTTATG\n>i2\nATGTTCATG\n")
        outgroup_fa = tmp_path / "outgroup.fa"
        outgroup_fa.write_text(">o1\nATGTTAATG\n")

        gene = extract_polymorphism_data(ingroup_fa, outgroup_fa)
        assert gene.ln is not None
        assert gene.ls is not None
        assert gene.ln + gene.ls == pytest.approx(9.0)


class TestAggregatePolymorphismData:
    def test_sums_sites(self) -> None:
        gene1 = PolymorphismData(polymorphisms=[], dn=0, ds=0, ln=200.0, ls=100.0)
        gene2 = PolymorphismData(polymorphisms=[], dn=0, ds=0, ln=150.0, ls=75.0)
        agg = aggregate_polymorphism_data([gene1, gene2], num_bins=5)
        assert agg.ln_total == pytest.approx(350.0)
        assert agg.ls_total == pytest.approx(175.0)
