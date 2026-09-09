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
def alpha_tg_ci_undefined() -> AlphaTGResult:
    """A single-gene result: point estimates are defined, CI fields are not."""
    return AlphaTGResult(
        alpha_tg=0.35,
        ni_tg=0.65,
        ci_low=None,
        ci_high=None,
        num_genes=1,
        dn_total=10,
        ds_total=20,
        pn_total=5,
        ps_total=15,
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
    "alpha_tg_ci_undefined",
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


def test_pretty_NA_appears_when_alpha_tg_ci_undefined(
    alpha_tg_ci_undefined: AlphaTGResult,
) -> None:
    """Targeted check: with one gene, AlphaTGResult.__str__ emits 'NA' for the CI."""
    pretty = str(alpha_tg_ci_undefined)
    assert "NA" in pretty, f"expected NA in pretty output; got:\n{pretty}"


def test_tsv_NA_appears_when_alpha_tg_ci_undefined(
    alpha_tg_ci_undefined: AlphaTGResult,
) -> None:
    """Targeted check: with one gene, AlphaTGResult TSV emits 'NA' for CI_low/CI_high."""
    tsv = format_result(alpha_tg_ci_undefined, OutputFormat.TSV)
    _, values = tsv.split("\n")
    assert "NA" in values.split("\t"), f"expected NA among TSV fields; got:\n{tsv}"


def test_batch_tsv_imputed(imputed_undefined: ImputedMKResult) -> None:
    """Batch TSV renders imputed results as one row per gene."""
    defined = ImputedMKResult(
        alpha=0.5,
        p_value=0.01,
        pn_neutral=3.0,
        pwd=1.0,
        dn=10,
        ds=5,
        pn_total=4,
        ps_total=8,
        cutoff=0.15,
        ci_method="bootstrap",
    )
    out = format_batch_results([("geneA", defined), ("geneB", imputed_undefined)], OutputFormat.TSV)
    lines = out.splitlines()
    assert lines[0].startswith("gene\tDn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha")
    assert len(lines) == 3
    assert lines[1].startswith("geneA\t10\t5\t4\t8\t1.00\t3.00\t0.500000\t0.01\t0.15\t")
    assert lines[2].startswith("geneB\t0\t0\t0\t0\t0.00\t0.00\tNA\t")


def test_batch_tsv_imputed_with_adjusted_pvalues(imputed_undefined: ImputedMKResult) -> None:
    """The adjusted p-value column, when present, sits directly after p_value."""
    defined = ImputedMKResult(
        alpha=0.5,
        p_value=0.01,
        pn_neutral=3.0,
        pwd=1.0,
        dn=10,
        ds=5,
        pn_total=4,
        ps_total=8,
        cutoff=0.15,
        ci_method="bootstrap",
    )
    out = format_batch_results(
        [("geneA", defined), ("geneB", imputed_undefined)],
        OutputFormat.TSV,
        adjusted_pvalues=[0.02, 0.9],
    )
    lines = out.splitlines()
    assert lines[0].startswith(
        "gene\tDn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tp_value_adjusted\tcutoff\t"
    )
    assert lines[1].startswith("geneA\t10\t5\t4\t8\t1.00\t3.00\t0.500000\t0.01\t0.02\t0.15\t")
    assert lines[2].startswith("geneB\t0\t0\t0\t0\t0.00\t0.00\tNA\t1\t0.9\t0.15\t")


def test_batch_json_polarized_adjusted_pvalue_nests_in_ingroup(
    polarized_undefined: PolarizedMKResult,
) -> None:
    """The adjusted p-value sits beside the p-value it adjusts, inside ``ingroup``."""
    out = format_batch_results(
        [("geneA", polarized_undefined), ("geneB", polarized_undefined)],
        OutputFormat.JSON,
        adjusted_pvalues=[0.02, 0.9],
    )
    data = json.loads(out)
    assert data["geneA"]["ingroup"]["p_value_adjusted"] == 0.02
    assert "p_value_adjusted" not in data["geneA"]
    assert data["geneB"]["ingroup"]["p_value_adjusted"] == 0.9


def test_batch_json_imputed_adjusted_pvalue_is_top_level(
    imputed_undefined: ImputedMKResult,
) -> None:
    """Unlike polarized, imputed's ``to_dict`` is flat, so the adjusted value stays top-level."""
    out = format_batch_results(
        [("geneA", imputed_undefined), ("geneB", imputed_undefined)],
        OutputFormat.JSON,
        adjusted_pvalues=[0.02, 0.9],
    )
    data = json.loads(out)
    assert data["geneA"]["p_value_adjusted"] == 0.02
    assert data["geneB"]["p_value_adjusted"] == 0.9


def test_batch_tsv_rejects_unknown_result_type(alpha_tg_undefined: AlphaTGResult) -> None:
    """A result type without a batch layout is an error, not a silent format change."""
    with pytest.raises(TypeError, match="Unknown result type"):
        format_batch_results([("geneA", alpha_tg_undefined)], OutputFormat.TSV)


@pytest.mark.parametrize(
    "first_fixture,second_fixture",
    [
        ("mk_undefined", "asymptotic_undefined"),
        ("asymptotic_undefined", "polarized_undefined"),
        ("polarized_undefined", "imputed_undefined"),
        ("imputed_undefined", "mk_undefined"),
    ],
)
def test_batch_tsv_rejects_mixed_result_types(
    first_fixture: str, second_fixture: str, request: pytest.FixtureRequest
) -> None:
    """A row whose type differs from the first row's is a caller error, not a dropped row."""
    first = request.getfixturevalue(first_fixture)
    second = request.getfixturevalue(second_fixture)
    with pytest.raises(TypeError, match="geneA.*geneB") as exc_info:
        format_batch_results([("geneA", first), ("geneB", second)], OutputFormat.TSV)
    message = str(exc_info.value)
    assert type(first).__name__ in message
    assert type(second).__name__ in message


def _defined_mk_result() -> MKResult:
    """A representative non-degenerate MKResult, reused wherever the exact counts don't matter."""
    return mk_test_from_counts(dn=10, ds=5, pn=4, ps=8)


def _mixed_result_types() -> list[tuple[str, MKResult | AsymptoticMKResult]]:
    """An MK/asymptotic/MK batch, none sharing the middle gene's type (issue #69's repro)."""
    return [
        ("geneA", _defined_mk_result()),
        (
            "geneB",
            AsymptoticMKResult(alpha_asymptotic=0.0, ci_low=0.0, ci_high=0.0, dn=1, ds=2),
        ),
        ("geneC", mk_test_from_counts(dn=6, ds=3, pn=2, ps=4)),
    ]


def test_batch_json_preserves_mixed_result_shapes() -> None:
    """JSON keys by gene name regardless of type, so a mixed batch is not rejected."""
    out = format_batch_results(_mixed_result_types(), OutputFormat.JSON)
    data = json.loads(out)
    assert set(data.keys()) == {"geneA", "geneB", "geneC"}


def test_batch_pretty_preserves_mixed_result_shapes() -> None:
    """Pretty output renders each result with its own str(), so a mixed batch is not rejected."""
    out = format_batch_results(_mixed_result_types(), OutputFormat.PRETTY)
    assert "=== geneA ===" in out
    assert "=== geneB ===" in out
    assert "=== geneC ===" in out


def _two_results_with_the_same_name() -> list[tuple[str, MKResult]]:
    return [
        ("geneA", _defined_mk_result()),
        ("geneA", mk_test_from_counts(dn=1, ds=1, pn=1, ps=1)),
    ]


def test_batch_json_rejects_duplicate_gene_names() -> None:
    """An object keyed by gene name cannot hold two results that share a name."""
    with pytest.raises(ValueError, match="geneA"):
        format_batch_results(_two_results_with_the_same_name(), OutputFormat.JSON)


@pytest.mark.parametrize("output_format", [OutputFormat.TSV, OutputFormat.PRETTY])
def test_batch_row_formats_keep_duplicate_gene_names(output_format: OutputFormat) -> None:
    """A row per result loses nothing to a repeated name, so it is not rejected."""
    out = format_batch_results(_two_results_with_the_same_name(), output_format)
    assert out.count("geneA") == 2
