"""Output formatters for MK test results."""

from __future__ import annotations

import json
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mkado.analysis.alpha_tg import AlphaTGResult
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult
    from mkado.analysis.mk_test import MKResult
    from mkado.analysis.polarized import PolarizedMKResult


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


def _omega_tsv_columns(result: MKResult | PolarizedMKResult) -> str:
    """Tab-separated ``Ln\\tLs\\tomega`` for omega-only result types (MK, polarized)."""
    return f"{_fmt_optional(result.ln)}\t{_fmt_optional(result.ls)}\t{_fmt_optional(result.omega)}"


_IMPUTED_TSV_HEADER = (
    "Dn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tcutoff\t"
    "Ln\tLs\tomega\tomega_a\tomega_na\t"
    "alpha_CI_low\talpha_CI_high\t"
    "omega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method"
)


def _imputed_tsv_row(result: ImputedMKResult) -> str:
    """One TSV row for an imputed result, in the column order of _IMPUTED_TSV_HEADER."""
    ci_method_str = result.ci_method if result.ci_method is not None else "NA"
    return (
        f"{result.dn}\t{result.ds}\t{result.pn_total}\t{result.ps_total}\t"
        f"{result.pwd:.2f}\t{result.pn_neutral:.2f}\t{_fmt_optional(result.alpha)}\t"
        f"{result.p_value:.6g}\t{result.cutoff}\t"
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
        return f"{_IMPUTED_TSV_HEADER}\n{_imputed_tsv_row(result)}"

    elif isinstance(result, MKResult):
        header = "Dn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega"
        ni_str = f"{result.ni:.6f}" if result.ni is not None else "NA"
        alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
        dos_str = f"{result.dos:.6f}" if result.dos is not None else "NA"
        values = (
            f"{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t{result.p_value:.6g}\t"
            f"{ni_str}\t{alpha_str}\t{dos_str}\t"
            f"{_omega_tsv_columns(result)}"
        )
        return f"{header}\n{values}"

    elif isinstance(result, PolarizedMKResult):
        header = "lineage\tDn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega"
        lines = [header]

        ni_str = f"{result.ni_ingroup:.6f}" if result.ni_ingroup is not None else "NA"
        alpha_str = f"{result.alpha_ingroup:.6f}" if result.alpha_ingroup is not None else "NA"
        dos_str = f"{result.dos_ingroup:.6f}" if result.dos_ingroup is not None else "NA"
        ln_s = _fmt_optional(result.ln)
        ls_s = _fmt_optional(result.ls)
        lines.append(
            f"ingroup\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
            f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
            f"{result.p_value_ingroup:.6g}\t{ni_str}\t{alpha_str}\t{dos_str}\t"
            f"{_omega_tsv_columns(result)}"
        )
        lines.append(
            f"outgroup\t{result.dn_outgroup}\t{result.ds_outgroup}\tNA\tNA\tNA\tNA\tNA\tNA\t"
            f"{ln_s}\t{ls_s}\tNA"
        )
        lines.append(
            f"unpolarized\t{result.dn_unpolarized}\t{result.ds_unpolarized}\t"
            f"NA\tNA\tNA\tNA\tNA\tNA\t{ln_s}\t{ls_s}\tNA"
        )
        return "\n".join(lines)

    elif isinstance(result, AsymptoticMKResult):
        ci_cols = (
            "\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high"
            "\tci_method\tsfs_mode"
        )
        ci_values = f"{result.ci_method}\t{result.sfs_mode}"
        if result.num_genes > 0:
            # Aggregated result
            header = (
                "Dn\tDs\tPn\tPs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tnum_genes\t"
                "Ln\tLs\tomega\tomega_a\tomega_na" + ci_cols
            )
            values = (
                f"{result.dn}\t{result.ds}\t{result.pn_total}\t{result.ps_total}\t"
                f"{result.alpha_asymptotic:.6f}\t{result.ci_low:.6f}\t{result.ci_high:.6f}\t"
                f"{result.model_type}\t{result.num_genes}\t"
                f"{_omega_full_tsv_columns(result)}\t"
                f"{_omega_a_na_ci_tsv_columns(result)}\t{ci_values}"
            )
        else:
            header = (
                "Dn\tDs\talpha_asymptotic\tCI_low\tCI_high\t"
                "Ln\tLs\tomega\tomega_a\tomega_na" + ci_cols
            )
            values = (
                f"{result.dn}\t{result.ds}\t{result.alpha_asymptotic:.6f}\t"
                f"{result.ci_low:.6f}\t{result.ci_high:.6f}\t"
                f"{_omega_full_tsv_columns(result)}\t"
                f"{_omega_a_na_ci_tsv_columns(result)}\t{ci_values}"
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
            f"{result.alpha_tg:.6f}\t{result.ni_tg:.6f}\t{result.ci_low:.6f}\t{result.ci_high:.6f}\t"
            f"{result.num_genes}\t"
            f"{_omega_full_tsv_columns(result)}\t"
            f"{_fmt_optional(result.omega_ci_low)}\t{_fmt_optional(result.omega_ci_high)}\t"
            f"{_omega_a_na_ci_tsv_columns(result)}\t{result.ci_method}"
        )
        return f"{header}\n{values}"

    else:
        raise TypeError(f"Unknown result type: {type(result)}")


def format_batch_results(
    results: list[tuple[str, MKResult | PolarizedMKResult | AsymptoticMKResult | ImputedMKResult]],
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
        data = {}
        for i, (name, result) in enumerate(results):
            if name in data:
                raise ValueError(
                    f"Duplicate gene name '{name}': JSON output is keyed by gene name. "
                    "Use tsv or pretty format."
                )
            result_dict = result.to_dict()
            if adjusted_pvalues is not None:
                result_dict["p_value_adjusted"] = adjusted_pvalues[i]
            data[name] = result_dict
        return json.dumps(data, indent=2)

    elif format == OutputFormat.TSV:
        from mkado.analysis.asymptotic import AsymptoticMKResult
        from mkado.analysis.imputed import ImputedMKResult
        from mkado.analysis.mk_test import MKResult
        from mkado.analysis.polarized import PolarizedMKResult

        if not results:
            return ""

        # Check type of first result
        _, first_result = results[0]
        if isinstance(first_result, MKResult):
            base_header = "gene\tDn\tDs\tPn\tPs\tp_value"
            if adjusted_pvalues is not None:
                base_header += "\tp_value_adjusted"
            header = base_header + "\tNI\talpha\tDoS\tLn\tLs\tomega"
            lines = [header]
            for i, (name, result) in enumerate(results):
                if isinstance(result, MKResult):
                    ni_str = f"{result.ni:.6f}" if result.ni is not None else "NA"
                    alpha_str = f"{result.alpha:.6f}" if result.alpha is not None else "NA"
                    dos_str = f"{result.dos:.6f}" if result.dos is not None else "NA"
                    base_row = (
                        f"{name}\t{result.dn}\t{result.ds}\t{result.pn}\t{result.ps}\t"
                        f"{result.p_value:.6g}"
                    )
                    if adjusted_pvalues is not None:
                        base_row += f"\t{adjusted_pvalues[i]:.6g}"
                    lines.append(
                        f"{base_row}\t{ni_str}\t{alpha_str}\t{dos_str}\t"
                        f"{_omega_tsv_columns(result)}"
                    )
            return "\n".join(lines)

        elif isinstance(first_result, AsymptoticMKResult):
            header = (
                "gene\tDn\tDs\talpha_asymptotic\tCI_low\tCI_high\tmodel\t"
                "Ln\tLs\tomega\tomega_a\tomega_na\t"
                "omega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\t"
                "ci_method\tsfs_mode"
            )
            lines = [header]
            for name, result in results:
                if isinstance(result, AsymptoticMKResult):
                    lines.append(
                        f"{name}\t{result.dn}\t{result.ds}\t"
                        f"{result.alpha_asymptotic:.6f}\t{result.ci_low:.6f}\t"
                        f"{result.ci_high:.6f}\t{result.model_type}\t"
                        f"{_omega_full_tsv_columns(result)}\t"
                        f"{_omega_a_na_ci_tsv_columns(result)}\t"
                        f"{result.ci_method}\t{result.sfs_mode}"
                    )
            return "\n".join(lines)

        elif isinstance(first_result, PolarizedMKResult):
            base_header = (
                "gene\tDn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\t"
                "Dn_outgroup\tDs_outgroup\tp_value"
            )
            if adjusted_pvalues is not None:
                base_header += "\tp_value_adjusted"
            header = base_header + "\tNI\talpha\tDoS\tLn\tLs\tomega"
            lines = [header]
            for i, (name, result) in enumerate(results):
                if isinstance(result, PolarizedMKResult):
                    ni_str = f"{result.ni_ingroup:.6f}" if result.ni_ingroup is not None else "NA"
                    alpha_str = (
                        f"{result.alpha_ingroup:.6f}" if result.alpha_ingroup is not None else "NA"
                    )
                    dos_str = (
                        f"{result.dos_ingroup:.6f}" if result.dos_ingroup is not None else "NA"
                    )
                    base_row = (
                        f"{name}\t{result.dn_ingroup}\t{result.ds_ingroup}\t"
                        f"{result.pn_ingroup}\t{result.ps_ingroup}\t"
                        f"{result.dn_outgroup}\t{result.ds_outgroup}\t"
                        f"{result.p_value_ingroup:.6g}"
                    )
                    if adjusted_pvalues is not None:
                        base_row += f"\t{adjusted_pvalues[i]:.6g}"
                    lines.append(
                        f"{base_row}\t{ni_str}\t{alpha_str}\t{dos_str}\t"
                        f"{_omega_tsv_columns(result)}"
                    )
            return "\n".join(lines)

        elif isinstance(first_result, ImputedMKResult):
            lines = ["gene\t" + _IMPUTED_TSV_HEADER]
            for name, result in results:
                if isinstance(result, ImputedMKResult):
                    lines.append(f"{name}\t{_imputed_tsv_row(result)}")
            return "\n".join(lines)

        else:
            raise TypeError(f"Unknown result type: {type(first_result)}")

    else:
        raise ValueError(f"Unknown format: {format}")
