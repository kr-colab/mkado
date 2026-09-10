"""Tests for the G-test and alpha confidence interval helpers."""

import math

import pytest
from scipy import stats

from mkado.analysis.statistics import alpha, confidence_interval_alpha, g_test


class TestGTest:
    """Tests for the G-test (log-likelihood ratio test)."""

    @pytest.mark.parametrize(
        "dn, ds, pn, ps",
        [
            (10, 5, 4, 8),
            (7, 17, 2, 42),
            (30, 10, 5, 40),
        ],
    )
    def test_matches_scipy_log_likelihood(self, dn: int, ds: int, pn: int, ps: int) -> None:
        """g_test matches scipy's log-likelihood chi2_contingency p-value."""
        scipy_p = stats.chi2_contingency(
            [[dn, ds], [pn, ps]], lambda_="log-likelihood", correction=False
        )[1]
        assert g_test(dn, ds, pn, ps) == pytest.approx(scipy_p, rel=1e-9)

    def test_equal_ratios_give_p_of_one(self) -> None:
        """Equal row/column ratios give G = 0 and p-value 1.0."""
        assert g_test(10, 10, 10, 10) == 1.0

    def test_all_zero_table_returns_one(self) -> None:
        """An all-zero table hits the total == 0 short-circuit and returns 1.0."""
        assert g_test(0, 0, 0, 0) == 1.0

    def test_zero_row_returns_one(self) -> None:
        """A table with an empty row has G = 0 and returns 1.0."""
        assert g_test(10, 5, 0, 0) == 1.0

    def test_zero_column_returns_one(self) -> None:
        """A table with an empty column has G = 0 and returns 1.0."""
        assert g_test(10, 0, 4, 0) == 1.0

    def test_strongly_significant_table_does_not_underflow(self) -> None:
        """A strongly significant table reports a nonzero p-value instead of underflowing."""
        dn, ds, pn, ps = 100, 10, 10, 100
        g_stat = stats.chi2_contingency(
            [[dn, ds], [pn, ps]], lambda_="log-likelihood", correction=False
        )[0]
        expected_p = stats.chi2.sf(g_stat, df=1)

        result = g_test(dn, ds, pn, ps)

        assert result > 0
        assert result == pytest.approx(expected_p, rel=1e-9)

    def test_returns_python_float(self) -> None:
        """g_test returns a plain Python float, not a numpy scalar."""
        assert type(g_test(10, 5, 4, 8)) is float


class TestConfidenceIntervalAlpha:
    """Tests for the delta-method confidence interval on alpha."""

    def test_hand_computed_interval(self) -> None:
        """The interval matches a hand-computed delta-method calculation."""
        r = 0.25
        var_r = r**2 * (1 / 10 + 1 / 5 + 1 / 4 + 1 / 8)
        se = math.sqrt(var_r)
        z = stats.norm.ppf(0.975)
        expected_lower = 0.75 - z * se
        expected_upper = 0.75 + z * se

        result = confidence_interval_alpha(10, 5, 4, 8)

        assert result is not None
        lower, upper = result
        assert lower == pytest.approx(expected_lower)
        assert upper == pytest.approx(expected_upper)
        assert lower == pytest.approx(0.3474313, abs=1e-6)
        assert upper == pytest.approx(1.1525687, abs=1e-6)

    @pytest.mark.parametrize(
        "dn, ds, pn, ps",
        [
            (10, 5, 4, 8),
            (6, 3, 2, 4),
            (20, 10, 5, 10),
        ],
    )
    def test_interval_contains_alpha(self, dn: int, ds: int, pn: int, ps: int) -> None:
        """The interval is centered on alpha and contains it."""
        a = alpha(dn, ds, pn, ps)
        result = confidence_interval_alpha(dn, ds, pn, ps)

        assert a is not None
        assert result is not None
        lower, upper = result
        assert lower < a < upper
        assert lower + upper == pytest.approx(2 * a)

    def test_higher_confidence_is_wider(self) -> None:
        """A higher confidence level widens the interval in both directions."""
        result_95 = confidence_interval_alpha(10, 5, 4, 8, confidence=0.95)
        result_99 = confidence_interval_alpha(10, 5, 4, 8, confidence=0.99)

        assert result_95 is not None
        assert result_99 is not None
        lower_95, upper_95 = result_95
        lower_99, upper_99 = result_99

        assert lower_99 < lower_95
        assert upper_99 > upper_95
        assert lower_95 == pytest.approx(0.3474313, abs=1e-6)
        assert upper_95 == pytest.approx(1.1525687, abs=1e-6)
        assert lower_99 == pytest.approx(0.2209351, abs=1e-6)
        assert upper_99 == pytest.approx(1.2790649, abs=1e-6)

    @pytest.mark.parametrize(
        "dn, ds, pn, ps",
        [
            (0, 5, 4, 8),
            (10, 5, 4, 0),
        ],
    )
    def test_undefined_alpha_returns_none(self, dn: int, ds: int, pn: int, ps: int) -> None:
        """When alpha itself is undefined, the interval is also None."""
        assert confidence_interval_alpha(dn, ds, pn, ps) is None

    @pytest.mark.parametrize(
        "dn, ds, pn, ps",
        [
            (10, 0, 4, 8),
            (10, 5, 0, 8),
        ],
    )
    def test_zero_ds_or_pn_uses_pseudo_counts(self, dn: int, ds: int, pn: int, ps: int) -> None:
        """A zero Ds or Pn still yields a finite interval via pseudo-counts."""
        a = alpha(dn, ds, pn, ps)
        result = confidence_interval_alpha(dn, ds, pn, ps)

        assert a == 1.0
        assert result is not None
        lower, upper = result
        assert math.isfinite(lower)
        assert math.isfinite(upper)
        assert lower < 1.0 < upper

    def test_default_confidence_is_95_percent(self) -> None:
        """The default confidence level matches an explicit 0.95."""
        assert confidence_interval_alpha(10, 5, 4, 8) == confidence_interval_alpha(
            10, 5, 4, 8, confidence=0.95
        )
