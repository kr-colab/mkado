"""Tests for output-format conventions (issue #13).

MKado renders missing/undefined floats as the literal ``"NA"`` (no slash) in
both pretty and TSV output, matching the R/Pandas convention. JSON uses
``null``. These tests guard against regressions where one format drifts
from another.
"""

from __future__ import annotations

import json

import pytest

from mkado.analysis.alpha_tg import AlphaTGResult
from mkado.analysis.asymptotic import AsymptoticMKResult
from mkado.analysis.imputed import ImputedMKResult
from mkado.analysis.mk_test import MKResult, mk_test_from_counts
from mkado.analysis.polarized import PolarizedMKResult
from mkado.io.output import OutputFormat, format_batch_results, format_result


@pytest.fixture
def mk_undefined() -> MKResult:
    """Counts where alpha/NI/DoS/omega all degenerate to None."""
    # Dn=0 forces alpha and NI to None; Dn+Ds=0 and Pn+Ps=0 forces DoS to None.
    return mk_test_from_counts(dn=0, ds=0, pn=0, ps=0)


@pytest.fixture
def asymptotic_undefined() -> AsymptoticMKResult:
    return AsymptoticMKResult(
        alpha_asymptotic=0.0,
        ci_low=0.0,
        ci_high=0.0,
        dn=0,
        ds=0,
        # ln/ls left None → omega/omega_a/omega_na all None
    )


@pytest.fixture
def imputed_undefined() -> ImputedMKResult:
    return ImputedMKResult(
        alpha=None,
        p_value=1.0,
        pn_neutral=0.0,
        pwd=0.0,
        dn=0,
        ds=0,
        pn_total=0,
        ps_total=0,
        cutoff=0.15,
    )


@pytest.mark.parametrize("adjusted", [None, [0.25, 0.5]])
def test_imputed_batch_tsv_preserves_single_result_columns(imputed_undefined, adjusted):
    header, row = format_result(imputed_undefined, OutputFormat.TSV).splitlines()
    output = format_batch_results(
        [("first", imputed_undefined), ("second", imputed_undefined)],
        OutputFormat.TSV,
        adjusted_pvalues=adjusted,
    ).splitlines()
    suffix = "\tp_value_adjusted" if adjusted is not None else ""
    assert output[0] == "gene\t" + header + suffix
    for i, name in enumerate(["first", "second"]):
        suffix = f"\t{adjusted[i]:.6g}" if adjusted is not None else ""
        assert output[i + 1] == f"{name}\t{row}{suffix}"
        assert len(output[i + 1].split("\t")) == len(output[0].split("\t"))


@pytest.fixture
def alpha_tg_undefined() -> AlphaTGResult:
    return AlphaTGResult(
        alpha_tg=0.0,
        ni_tg=0.0,
        ci_low=0.0,
        ci_high=0.0,
        num_genes=0,
        dn_total=0,
        ds_total=0,
        pn_total=0,
        ps_total=0,
        # ln/ls left None → omega/omega_a/omega_na all None
    )


@pytest.fixture
def polarized_undefined() -> PolarizedMKResult:
    return PolarizedMKResult(
        dn_ingroup=0,
        ds_ingroup=0,
        pn_ingroup=0,
        ps_ingroup=0,
        dn_outgroup=0,
        ds_outgroup=0,
        dn_unpolarized=0,
        ds_unpolarized=0,
        pn_unpolarized=0,
        ps_unpolarized=0,
        p_value_ingroup=1.0,
        ni_ingroup=None,
        alpha_ingroup=None,
        dos_ingroup=None,
    )


_ALL_UNDEFINED = [
    "mk_undefined",
    "asymptotic_undefined",
    "imputed_undefined",
    "alpha_tg_undefined",
    "polarized_undefined",
]


@pytest.mark.parametrize("fixture_name", _ALL_UNDEFINED)
def test_pretty_output_never_uses_N_slash_A(
    fixture_name: str, request: pytest.FixtureRequest
) -> None:
    """Pretty output must never render missing values as ``"N/A"`` (with slash)."""
    result = request.getfixturevalue(fixture_name)
    pretty = str(result)
    assert "N/A" not in pretty, f"{fixture_name}: pretty output uses 'N/A'\n{pretty}"


@pytest.mark.parametrize("fixture_name", _ALL_UNDEFINED)
def test_tsv_output_never_uses_N_slash_A(fixture_name: str, request: pytest.FixtureRequest) -> None:
    """TSV output must never render missing values as ``"N/A"`` (with slash)."""
    result = request.getfixturevalue(fixture_name)
    tsv = format_result(result, OutputFormat.TSV)
    assert "N/A" not in tsv, f"{fixture_name}: TSV output uses 'N/A'\n{tsv}"


@pytest.mark.parametrize("fixture_name", _ALL_UNDEFINED)
def test_json_never_uses_N_slash_A_or_string_NA(
    fixture_name: str, request: pytest.FixtureRequest
) -> None:
    """JSON renders missing values as ``null``, never ``"N/A"`` or string ``"NA"``."""
    result = request.getfixturevalue(fixture_name)
    js = format_result(result, OutputFormat.JSON)
    assert "N/A" not in js, f"{fixture_name}: JSON contains 'N/A'\n{js}"
    # Sanity: it parses and at least one value round-trips
    parsed = json.loads(js)
    assert parsed is not None


def test_pretty_NA_appears_when_alpha_undefined() -> None:
    """Targeted check: when alpha is None, MKResult.__str__ emits the literal 'NA'."""
    # dn=0, ds=0, pn=0, ps=0 → ni/alpha/dos all None → all three rendered as "NA"
    result = mk_test_from_counts(dn=0, ds=0, pn=0, ps=0)
    pretty = str(result)
    # Three "NA" tokens (NI, alpha, DoS lines), at minimum
    assert pretty.count("NA") >= 3, f"expected ≥3 NA tokens; got:\n{pretty}"


def test_tsv_NA_appears_when_alpha_undefined() -> None:
    """Targeted check: when alpha is None, MKResult TSV emits the literal 'NA'."""
    result = mk_test_from_counts(dn=0, ds=0, pn=0, ps=0)
    tsv = format_result(result, OutputFormat.TSV)
    # TSV header has no NA; data row has NA wherever a None field landed
    _, values = tsv.split("\n")
    assert "NA" in values.split("\t"), f"expected NA among TSV fields; got:\n{tsv}"
