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


def test_batch_tsv_rejects_unknown_result_type(alpha_tg_undefined: AlphaTGResult) -> None:
    """A result type without a batch layout is an error, not a silent format change."""
    with pytest.raises(TypeError, match="Unknown result type"):
        format_batch_results([("geneA", alpha_tg_undefined)], OutputFormat.TSV)


def _two_results_with_the_same_name() -> list[tuple[str, MKResult]]:
    return [
        ("geneA", mk_test_from_counts(dn=10, ds=5, pn=4, ps=8)),
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


@pytest.fixture
def two_mk_results() -> list[tuple[str, MKResult]]:
    """Two distinct records make row loss and name mixups visible."""
    return [
        ("geneA", mk_test_from_counts(dn=10, ds=5, pn=4, ps=8)),
        ("geneB", mk_test_from_counts(dn=8, ds=2, pn=3, ps=9)),
    ]


@pytest.mark.parametrize("output_format", list(OutputFormat))
def test_batch_formats_keep_both_gene_names(
    two_mk_results: list[tuple[str, MKResult]], output_format: OutputFormat
) -> None:
    output = format_batch_results(two_mk_results, output_format)
    assert "geneA" in output
    assert "geneB" in output
    if output_format == OutputFormat.TSV:
        assert len(output.splitlines()) == 3
    elif output_format == OutputFormat.JSON:
        assert list(json.loads(output)) == ["geneA", "geneB"]
    else:
        assert output.count("=== ") == 2


@pytest.mark.parametrize("output_format", list(OutputFormat))
def test_batch_adjusted_pvalues_are_attached_to_matching_results(
    two_mk_results: list[tuple[str, MKResult]], output_format: OutputFormat
) -> None:
    output = format_batch_results(two_mk_results, output_format, adjusted_pvalues=[0.02, 0.1])
    if output_format == OutputFormat.JSON:
        parsed = json.loads(output)
        assert parsed["geneA"]["p_value_adjusted"] == 0.02
        assert parsed["geneB"]["p_value_adjusted"] == 0.1
    elif output_format == OutputFormat.TSV:
        rows = [line.split("\t") for line in output.splitlines()]
        adjusted_index = rows[0].index("p_value_adjusted")
        assert rows[1][adjusted_index] == "0.02"
        assert rows[2][adjusted_index] == "0.1"
    else:
        blocks = {block.split(" ===", 1)[0]: block for block in output.split("=== ")[1:]}
        assert "p-value (BH adj):     0.02" in blocks["geneA"]
        assert "p-value (BH adj):     0.1" in blocks["geneB"]


@pytest.mark.parametrize("output_format", list(OutputFormat))
@pytest.mark.parametrize("adjusted_pvalues", [[], [0.01], [0.01, 0.02, 0.03]])
def test_batch_rejects_adjusted_pvalue_length_mismatch(
    two_mk_results: list[tuple[str, MKResult]],
    adjusted_pvalues: list[float],
    output_format: OutputFormat,
) -> None:
    with pytest.raises(ValueError, match="same length"):
        format_batch_results(
            two_mk_results, output_format, adjusted_pvalues=adjusted_pvalues
        )


@pytest.mark.parametrize(
    ("output_format", "expected"),
    [(OutputFormat.PRETTY, ""), (OutputFormat.TSV, ""), (OutputFormat.JSON, "{}")],
)
def test_batch_empty_results(output_format: OutputFormat, expected: str) -> None:
    assert format_batch_results([], output_format) == expected


def test_batch_json_preserves_mixed_result_shapes(
    mk_undefined: MKResult,
    polarized_undefined: PolarizedMKResult,
    asymptotic_undefined: AsymptoticMKResult,
) -> None:
    output = format_batch_results(
        [
            ("standard", mk_undefined),
            ("polarized", polarized_undefined),
            ("asymptotic", asymptotic_undefined),
        ],
        OutputFormat.JSON,
    )
    parsed = json.loads(output)
    assert parsed["standard"]["dn"] == 0
    assert parsed["polarized"]["ingroup"]["dn"] == 0
    assert parsed["asymptotic"]["alpha_asymptotic"] == 0.0
