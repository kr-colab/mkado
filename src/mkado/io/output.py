"""Output formatters for MK test results."""

from __future__ import annotations

import json
from dataclasses import dataclass
from enum import Enum
from functools import cache
from typing import TYPE_CHECKING, Any, Callable, NoReturn

if TYPE_CHECKING:
    from mkado.analysis.alpha_tg import AlphaTGResult
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult

    # The four result types that carry a per-gene batch layout (TSV/JSON/pretty).
    # AlphaTGResult has no batch layout, so it is deliberately excluded here.
    BatchResult = MKResult | PolarizedMKResult | AsymptoticMKResult | ImputedMKResult


class OutputFormat(Enum):
    """Supported output formats."""

    PRETTY = "pretty"
    TSV = "tsv"
    JSON = "json"


def _fmt_optional(value: float | None, precision: int = 6) -> str:
    """Format optional float as ``"NA"`` if ``None``, else ``f"{value:.{precision}f}"``."""
    if value is None:
        return "NA"
    return f"{value:.{precision}f}"


def _omega_full_tsv_columns(
    result: AsymptoticMKResult | ImputedMKResult | AlphaTGResult,
) -> str:
    """Tab-separated ``Ln\\tLs\\tomega\\tomega_a\\tomega_na`` for full-decomposition results."""
    return (
        f"{_fmt_optional(result.ln)}\t{_fmt_optional(result.ls)}\t"
        f"{_fmt_optional(result.omega)}\t{_fmt_optional(result.omega_a)}\t"
        f"{_fmt_optional(result.omega_na)}"
    )


def _omega_a_na_ci_tsv_columns(
    result: AsymptoticMKResult | AlphaTGResult | ImputedMKResult,
) -> str:
    """Tab-separated CI columns for omega_a and omega_na."""
    return (
        f"{_fmt_optional(result.omega_a_ci_low)}\t{_fmt_optional(result.omega_a_ci_high)}\t"
        f"{_fmt_optional(result.omega_na_ci_low)}\t{_fmt_optional(result.omega_na_ci_high)}"
    )


_INGROUP_STATS_HEADER = "NI\talpha\tDoS\tLn\tLs\tomega"
"""The statistics columns of a standard or polarized row, after the p-value columns."""


def _ingroup_stats_tsv_columns(
    ni: float | None, alpha: float | None, dos: float | None, result: MKResult | PolarizedMKResult
) -> str:
    """The values of _INGROUP_STATS_HEADER for the standard and polarized tests."""
    return (
        f"{_fmt_optional(ni)}\t{_fmt_optional(alpha)}\t{_fmt_optional(dos)}\t"
        f"{_fmt_optional(result.ln)}\t{_fmt_optional(result.ls)}\t{_fmt_optional(result.omega)}"
    )


def _with_adjusted_header(header: str, include_adjusted: bool) -> str:
    """``header`` followed by the adjusted p-value column, where one is asked for."""
    return f"{header}\tp_value_adjusted" if include_adjusted else header


def _with_adjusted(columns: str, adjusted_p_value: float | None) -> str:
    """``columns`` followed by the adjusted p-value column, where one is given."""
    if adjusted_p_value is None:
        return columns
    return f"{columns}\t{adjusted_p_value:.6g}"


def _mk_tsv_header(include_adjusted: bool = False) -> str:
    """Column header for a standard MK row, matching _mk_tsv_row's order."""
    header = _with_adjusted_header("Dn\tDs\tPn\tPs\tp_value", include_adjusted)
    return f"{header}\t{_INGROUP_STATS_HEADER}"


def _mk_tsv_row(result: MKResult, adjusted_p_value: float | None = None) -> str:
    """One TSV row for a standard MK result, in the column order of _mk_tsv_header()."""
    counts = f"{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t{result.p_value:.6g}"
    stats = _ingroup_stats_tsv_columns(result.ni, result.alpha, result.dos, result)
    return f"{_with_adjusted(counts, adjusted_p_value)}\t{stats}"


def _polarized_batch_tsv_header(include_adjusted: bool = False) -> str:
    """Column header for a polarized batch row: one row per gene, both lineages' counts."""
    header = _with_adjusted_header(
        "Dn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\tDn_outgroup\tDs_outgroup\tp_value",
        include_adjusted,
    )
    return f"{header}\t{_INGROUP_STATS_HEADER}"


def _polarized_batch_tsv_row(
    result: PolarizedMKResult, adjusted_p_value: float | None = None
) -> str:
    """One polarized batch row, in the column order of _polarized_batch_tsv_header()."""
    counts = (
        f"{result.dn_ingroup}\t{result.ds_ingroup}\t{result.pn_ingroup}\t{result.ps_ingroup}\t"
        f"{result.dn_outgroup}\t{result.ds_outgroup}\t{result.p_value_ingroup:.6g}"
    )
    stats = _ingroup_stats_tsv_columns(
        result.ni_ingroup, result.alpha_ingroup, result.dos_ingroup, result
    )
    return f"{_with_adjusted(counts, adjusted_p_value)}\t{stats}"


_ASYMPTOTIC_OMEGA_HEADER = (
    "Ln\tLs\tomega\tomega_a\tomega_na\t"
    "omega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\tsfs_mode"
)
"""The omega columns, their CIs, ``ci_method`` and ``sfs_mode`` of an asymptotic row."""


def _asymptotic_omega_tsv_columns(result: AsymptoticMKResult) -> str:
    """The values of _ASYMPTOTIC_OMEGA_HEADER for an asymptotic result."""
    return (
        f"{_omega_full_tsv_columns(result)}\t{_omega_a_na_ci_tsv_columns(result)}\t"
        f"{result.ci_method}\t{result.sfs_mode}"
    )


def _asymptotic_tsv_header(include_adjusted: bool = False) -> str:
    """Column header for a per-gene asymptotic row. The test has no p-value, so
    ``include_adjusted`` adds no column."""
    return f"Dn\tDs\talpha_asymptotic\tCI_low\tCI_high\tmodel\t{_ASYMPTOTIC_OMEGA_HEADER}"


def _asymptotic_tsv_row(result: AsymptoticMKResult, adjusted_p_value: float | None = None) -> str:
    """One per-gene asymptotic row, in the column order of _asymptotic_tsv_header()."""
    return (
        f"{result.dn}\t{result.ds}\t{result.alpha_asymptotic:.6f}\t{result.ci_low:.6f}\t"
        f"{result.ci_high:.6f}\t{result.model_type}\t{_asymptotic_omega_tsv_columns(result)}"
    )


def _imputed_tsv_header(include_adjusted: bool = False) -> str:
    """Column header for an imputed result row, matching _imputed_tsv_row's order."""
    header = _with_adjusted_header(
        "Dn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value", include_adjusted
    )
    return (
        f"{header}\tcutoff\t"
        "Ln\tLs\tomega\tomega_a\tomega_na\t"
        "alpha_CI_low\talpha_CI_high\t"
        "omega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method"
    )


def _imputed_tsv_row(result: ImputedMKResult, adjusted_p_value: float | None = None) -> str:
    """One TSV row for an imputed result, in the column order of _imputed_tsv_header()."""
    ci_method_str = result.ci_method if result.ci_method is not None else "NA"
    row = (
        f"{result.dn}\t{result.ds}\t{result.pn_total}\t{result.ps_total}\t"
        f"{result.pwd:.2f}\t{result.pn_neutral:.2f}\t{_fmt_optional(result.alpha)}\t"
        f"{result.p_value:.6g}"
    )
    return (
        f"{_with_adjusted(row, adjusted_p_value)}\t{result.cutoff}\t"
        f"{_omega_full_tsv_columns(result)}\t"
        f"{_fmt_optional(result.alpha_ci_low)}\t{_fmt_optional(result.alpha_ci_high)}\t"
        f"{_omega_a_na_ci_tsv_columns(result)}\t"
        f"{ci_method_str}"
    )


def format_result(
    result: MKResult | PolarizedMKResult | AsymptoticMKResult | AlphaTGResult | ImputedMKResult,
    format: OutputFormat = OutputFormat.PRETTY,
) -> str:
    """Format MK test results for output.

    Args:
        result: MK test result object
        format: Output format

    Returns:
        Formatted string representation

    Raises:
        ValueError: If the requested output format is unknown.
        TypeError: If TSV output does not support the result type.
    """
    if format == OutputFormat.PRETTY:
        return str(result)

    elif format == OutputFormat.JSON:
        return json.dumps(result.to_dict(), indent=2)

    elif format == OutputFormat.TSV:
        return _format_tsv(result)

    else:
        raise ValueError(f"Unknown format: {format}")


def _format_tsv(
    result: MKResult | PolarizedMKResult | AsymptoticMKResult | AlphaTGResult | ImputedMKResult,
) -> str:
    """Format results as tab-separated values."""
    from mkado.analysis.alpha_tg import AlphaTGResult
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult

    if isinstance(result, ImputedMKResult):
        return f"{_imputed_tsv_header()}\n{_imputed_tsv_row(result)}"

    elif isinstance(result, MKResult):
        return f"{_mk_tsv_header()}\n{_mk_tsv_row(result)}"

    elif isinstance(result, PolarizedMKResult):
        # One row per lineage; only the ingroup lineage has polymorphism and statistics.
        ln_ls = f"{_fmt_optional(result.ln)}\t{_fmt_optional(result.ls)}"
        no_stats = "NA\tNA\tNA\tNA\tNA\tNA"
        stats = _ingroup_stats_tsv_columns(
            result.ni_ingroup, result.alpha_ingroup, result.dos_ingroup, result
        )
        return "\n".join(
            [
                f"lineage\t{_mk_tsv_header()}",
                f"ingroup\t{result.dn_ingroup}\t{result.ds_ingroup}\t{result.pn_ingroup}\t"
                f"{result.ps_ingroup}\t{result.p_value_ingroup:.6g}\t{stats}",
                f"outgroup\t{result.dn_outgroup}\t{result.ds_outgroup}\t{no_stats}\t{ln_ls}\tNA",
                f"unpolarized\t{result.dn_unpolarized}\t{result.ds_unpolarized}\t{no_stats}\t"
                f"{ln_ls}\tNA",
            ]
        )

    elif isinstance(result, AsymptoticMKResult):
        if result.num_genes == 0:
            return f"{_asymptotic_tsv_header()}\n{_asymptotic_tsv_row(result)}"
        # An aggregated result also reports the pooled polymorphism counts and gene count.
        header = (
            "Dn\tDs\tPn\tPs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tnum_genes\t"
            + _ASYMPTOTIC_OMEGA_HEADER
        )
        values = (
            f"{result.dn}\t{result.ds}\t{result.pn_total}\t{result.ps_total}\t"
            f"{result.alpha_asymptotic:.6f}\t{result.ci_low:.6f}\t{result.ci_high:.6f}\t"
            f"{result.model_type}\t{result.num_genes}\t{_asymptotic_omega_tsv_columns(result)}"
        )
        return f"{header}\n{values}"

    elif isinstance(result, AlphaTGResult):
        header = (
            "Dn\tDs\tPn\tPs\talpha_TG\tNI_TG\tCI_low\tCI_high\tnum_genes\t"
            "Ln\tLs\tomega\tomega_a\tomega_na\t"
            "omega_CI_low\tomega_CI_high\t"
            "omega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method"
        )
        values = (
            f"{result.dn_total}\t{result.ds_total}\t{result.pn_total}\t{result.ps_total}\t"
            f"{result.alpha_tg:.6f}\t{result.ni_tg:.6f}\t"
            f"{_fmt_optional(result.ci_low)}\t{_fmt_optional(result.ci_high)}\t"
            f"{result.num_genes}\t"
            f"{_omega_full_tsv_columns(result)}\t"
            f"{_fmt_optional(result.omega_ci_low)}\t{_fmt_optional(result.omega_ci_high)}\t"
            f"{_omega_a_na_ci_tsv_columns(result)}\t{result.ci_method}"
        )
        return f"{header}\n{values}"

    else:
        raise TypeError(f"Unknown result type: {type(result)}")


@dataclass(frozen=True)
class _BatchLayout:
    """How one result type is written as a batch TSV row, and the p-value that the
    Benjamini-Hochberg adjustment uses (None for a type with no p-value)."""

    header: Callable[[bool], str]
    row: Callable[[Any, float | None], str]
    p_value: Callable[[Any], float | None]


@cache
def _batch_layouts() -> dict[type, _BatchLayout]:
    """The batch layout of each result type that has one. AlphaTGResult has none."""
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult

    return {
        MKResult: _BatchLayout(_mk_tsv_header, _mk_tsv_row, lambda r: r.p_value),
        PolarizedMKResult: _BatchLayout(
            _polarized_batch_tsv_header, _polarized_batch_tsv_row, lambda r: r.p_value_ingroup
        ),
        AsymptoticMKResult: _BatchLayout(
            _asymptotic_tsv_header, _asymptotic_tsv_row, lambda r: None
        ),
        ImputedMKResult: _BatchLayout(_imputed_tsv_header, _imputed_tsv_row, lambda r: r.p_value),
    }


def _batch_layout(result: object) -> _BatchLayout:
    """The batch layout of ``result``'s type; TypeError for a type with none."""
    layout = _batch_layouts().get(type(result))
    if layout is None:
        raise TypeError(f"Unknown result type: {type(result)}")
    return layout


def batch_p_value(result: BatchResult) -> float | None:
    """The p-value that the batch Benjamini-Hochberg adjustment uses, or None where the
    result type has no p-value.

    Raises:
        TypeError: If the result type has no batch layout.
    """
    return _batch_layout(result).p_value(result)


def _raise_mixed_batch_type(
    first_name: str, first_type_name: str, name: str, type_name: str
) -> NoReturn:
    """Fail loudly on a batch TSV row whose type differs from the first row's.

    A single TSV table has one column layout for every row; a mismatched
    type would otherwise be silently skipped, dropping the gene with no
    error and no trace in the output.
    """
    raise TypeError(
        f"Batch TSV requires one result type per file: gene '{first_name}' is "
        f"{first_type_name}, but gene '{name}' is {type_name}. "
        "Use JSON output for a batch that mixes result types."
    )


def format_batch_results(
    results: list[tuple[str, BatchResult]],
    format: OutputFormat = OutputFormat.PRETTY,
    adjusted_pvalues: list[float] | None = None,
) -> str:
    """Format multiple MK test results for batch output.

    Args:
        results: List of (name, result) tuples
        format: Output format
        adjusted_pvalues: Optional list of Benjamini-Hochberg adjusted p-values
            (must be same length as results if provided)

    Returns:
        Formatted string representation

    Raises:
        ValueError: If ``adjusted_pvalues`` does not have one value per result,
            if JSON output meets a duplicate gene name, or if the requested
            output format is unknown.
        TypeError: If TSV output meets a result type with no batch layout, or
            a batch that mixes result types.
    """
    if adjusted_pvalues is not None and len(adjusted_pvalues) != len(results):
        raise ValueError("adjusted_pvalues must have the same length as results")

    if format == OutputFormat.PRETTY:
        lines = []
        for i, (name, result) in enumerate(results):
            lines.append(f"=== {name} ===")
            lines.append(str(result))
            if adjusted_pvalues is not None:
                lines.append(f"  p-value (BH adj):     {adjusted_pvalues[i]:.6g}")
            lines.append("")
        return "\n".join(lines)

    elif format == OutputFormat.JSON:
        from mkado.analysis.polarized import PolarizedMKResult

        data = {}
        for i, (name, result) in enumerate(results):
            if name in data:
                raise ValueError(
                    f"Duplicate gene name '{name}': JSON output is keyed by gene name. "
                    "Use tsv or pretty format."
                )
            result_dict = result.to_dict()
            if adjusted_pvalues is not None:
                if isinstance(result, PolarizedMKResult):
                    # Nest next to the p-value it adjusts (ingroup.p_value) rather than
                    # sitting as an ambiguous sibling of ingroup/outgroup/unpolarized.
                    result_dict["ingroup"]["p_value_adjusted"] = adjusted_pvalues[i]
                else:
                    result_dict["p_value_adjusted"] = adjusted_pvalues[i]
            data[name] = result_dict
        return json.dumps(data, indent=2)

    elif format == OutputFormat.TSV:
        if not results:
            return ""

        # Every row must share the first row's type: a TSV table has one column
        # layout for all its rows, so a mismatch is a caller error, not a row to skip.
        first_name, first_result = results[0]
        layout = _batch_layout(first_result)
        for name, result in results[1:]:
            if type(result) is not type(first_result):
                _raise_mixed_batch_type(
                    first_name, type(first_result).__name__, name, type(result).__name__
                )

        lines = [f"gene\t{layout.header(adjusted_pvalues is not None)}"]
        for i, (name, result) in enumerate(results):
            adjusted = adjusted_pvalues[i] if adjusted_pvalues is not None else None
            lines.append(f"{name}\t{layout.row(result, adjusted)}")
        return "\n".join(lines)

    else:
        raise ValueError(f"Unknown format: {format}")
