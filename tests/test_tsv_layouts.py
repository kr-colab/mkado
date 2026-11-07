"""Exact TSV text for each result type, in single and batch output.

These tests pin the column layouts that the single and batch formatters share, so a change to
a shared layout shows here as a changed string.
"""

from __future__ import annotations

import pytest

from mkado.analysis.alpha_tg import AlphaTGResult
from mkado.analysis.asymptotic import AsymptoticMKResult
from mkado.analysis.imputed import ImputedMKResult
from mkado.analysis.mk_test import mk_test_from_counts
from mkado.analysis.polarized import PolarizedMKResult
from mkado.io.output import OutputFormat, batch_p_value, format_batch_results, format_result

TSV = OutputFormat.TSV

MK = mk_test_from_counts(dn=6, ds=8, pn=1, ps=8, ln=573.0, ls=195.0)
POLARIZED = PolarizedMKResult(
    dn_ingroup=5,
    ds_ingroup=7,
    pn_ingroup=2,
    ps_ingroup=9,
    dn_outgroup=3,
    ds_outgroup=4,
    dn_unpolarized=1,
    ds_unpolarized=2,
    pn_unpolarized=0,
    ps_unpolarized=1,
    p_value_ingroup=0.04321,
    ni_ingroup=0.3889,
    alpha_ingroup=0.6111,
    dos_ingroup=0.2348,
    ln=500.0,
    ls=150.0,
    omega=0.2143,
)
ASYMPTOTIC_PER_GENE = AsymptoticMKResult(
    alpha_asymptotic=0.4567,
    ci_low=0.1234,
    ci_high=0.789,
    dn=20,
    ds=30,
    model_type="linear",
    ln=900.0,
    ls=300.0,
    omega=0.2222,
    omega_a=0.1015,
    omega_na=0.1207,
    ci_method="bootstrap",
    sfs_mode="xi",
)
ASYMPTOTIC_AGGREGATED = AsymptoticMKResult(
    alpha_asymptotic=0.5,
    ci_low=0.4,
    ci_high=0.6,
    dn=200,
    ds=300,
    num_genes=12,
    pn_total=150,
    ps_total=400,
    model_type="exponential",
    ln=9000.0,
    ls=3000.0,
    omega=0.2222,
    omega_a=0.1111,
    omega_na=0.1111,
)
IMPUTED = ImputedMKResult(
    alpha=0.25,
    p_value=0.0123,
    pn_neutral=7.5,
    pwd=2.5,
    dn=12,
    ds=20,
    pn_total=10,
    ps_total=30,
    cutoff=0.15,
)
ALPHA_TG = AlphaTGResult(
    alpha_tg=0.35,
    ni_tg=0.65,
    ci_low=0.1,
    ci_high=0.5,
    num_genes=4,
    dn_total=10,
    ds_total=20,
    pn_total=5,
    ps_total=15,
)
ADJUSTED = [0.0456, 0.0789]

SINGLE_MK = "Dn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega\n6\t8\t1\t8\t0.175957\t0.166667\t0.833333\t0.317460\t573.000000\t195.000000\t0.255236"
SINGLE_POLARIZED = "lineage\tDn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega\ningroup\t5\t7\t2\t9\t0.04321\t0.388900\t0.611100\t0.234800\t500.000000\t150.000000\t0.214300\noutgroup\t3\t4\tNA\tNA\tNA\tNA\tNA\tNA\t500.000000\t150.000000\tNA\nunpolarized\t1\t2\tNA\tNA\tNA\tNA\tNA\tNA\t500.000000\t150.000000\tNA"
SINGLE_ASYMPTOTIC_AGGREGATED = "Dn\tDs\tPn\tPs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tnum_genes\tLn\tLs\tomega\tomega_a\tomega_na\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\tsfs_mode\n200\t300\t150\t400\t0.500000\t0.400000\t0.600000\texponential\t12\t9000.000000\t3000.000000\t0.222200\t0.111100\t0.111100\t0.088880\t0.133320\t0.088880\t0.133320\tmonte-carlo\tat"
SINGLE_IMPUTED = "Dn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tcutoff\tLn\tLs\tomega\tomega_a\tomega_na\talpha_CI_low\talpha_CI_high\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\n12\t20\t10\t30\t2.50\t7.50\t0.250000\t0.0123\t0.15\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
SINGLE_ALPHA_TG = "Dn\tDs\tPn\tPs\talpha_TG\tNI_TG\tCI_low\tCI_high\tnum_genes\tLn\tLs\tomega\tomega_a\tomega_na\tomega_CI_low\tomega_CI_high\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\n10\t20\t5\t15\t0.350000\t0.650000\t0.100000\t0.500000\t4\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tbootstrap"
BATCH_MK = "gene\tDn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega\ng1\t6\t8\t1\t8\t0.175957\t0.166667\t0.833333\t0.317460\t573.000000\t195.000000\t0.255236\ng2\t6\t8\t1\t8\t0.175957\t0.166667\t0.833333\t0.317460\t573.000000\t195.000000\t0.255236"
BATCH_MK_ADJUSTED = "gene\tDn\tDs\tPn\tPs\tp_value\tp_value_adjusted\tNI\talpha\tDoS\tLn\tLs\tomega\ng1\t6\t8\t1\t8\t0.175957\t0.0456\t0.166667\t0.833333\t0.317460\t573.000000\t195.000000\t0.255236\ng2\t6\t8\t1\t8\t0.175957\t0.0789\t0.166667\t0.833333\t0.317460\t573.000000\t195.000000\t0.255236"
BATCH_POLARIZED = "gene\tDn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\tDn_outgroup\tDs_outgroup\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega\ng1\t5\t7\t2\t9\t3\t4\t0.04321\t0.388900\t0.611100\t0.234800\t500.000000\t150.000000\t0.214300\ng2\t5\t7\t2\t9\t3\t4\t0.04321\t0.388900\t0.611100\t0.234800\t500.000000\t150.000000\t0.214300"
BATCH_POLARIZED_ADJUSTED = "gene\tDn_ingroup\tDs_ingroup\tPn_ingroup\tPs_ingroup\tDn_outgroup\tDs_outgroup\tp_value\tp_value_adjusted\tNI\talpha\tDoS\tLn\tLs\tomega\ng1\t5\t7\t2\t9\t3\t4\t0.04321\t0.0456\t0.388900\t0.611100\t0.234800\t500.000000\t150.000000\t0.214300\ng2\t5\t7\t2\t9\t3\t4\t0.04321\t0.0789\t0.388900\t0.611100\t0.234800\t500.000000\t150.000000\t0.214300"
BATCH_ASYMPTOTIC = "gene\tDn\tDs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tLn\tLs\tomega\tomega_a\tomega_na\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\tsfs_mode\ng1\t20\t30\t0.456700\t0.123400\t0.789000\tlinear\t900.000000\t300.000000\t0.222200\t0.101500\t0.120700\t0.027419\t0.175316\t0.046884\t0.194781\tbootstrap\txi\ng2\t20\t30\t0.456700\t0.123400\t0.789000\tlinear\t900.000000\t300.000000\t0.222200\t0.101500\t0.120700\t0.027419\t0.175316\t0.046884\t0.194781\tbootstrap\txi"
BATCH_IMPUTED = "gene\tDn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tcutoff\tLn\tLs\tomega\tomega_a\tomega_na\talpha_CI_low\talpha_CI_high\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\ng1\t12\t20\t10\t30\t2.50\t7.50\t0.250000\t0.0123\t0.15\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\ng2\t12\t20\t10\t30\t2.50\t7.50\t0.250000\t0.0123\t0.15\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
BATCH_IMPUTED_ADJUSTED = "gene\tDn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tp_value_adjusted\tcutoff\tLn\tLs\tomega\tomega_a\tomega_na\talpha_CI_low\talpha_CI_high\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\ng1\t12\t20\t10\t30\t2.50\t7.50\t0.250000\t0.0123\t0.0456\t0.15\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\ng2\t12\t20\t10\t30\t2.50\t7.50\t0.250000\t0.0123\t0.0789\t0.15\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
SINGLE_ASYMPTOTIC_PER_GENE = "Dn\tDs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tLn\tLs\tomega\tomega_a\tomega_na\tomega_a_CI_low\tomega_a_CI_high\tomega_na_CI_low\tomega_na_CI_high\tci_method\tsfs_mode\n20\t30\t0.456700\t0.123400\t0.789000\tlinear\t900.000000\t300.000000\t0.222200\t0.101500\t0.120700\t0.027419\t0.175316\t0.046884\t0.194781\tbootstrap\txi"


@pytest.mark.parametrize(
    ("result", "expected"),
    [
        (MK, SINGLE_MK),
        (POLARIZED, SINGLE_POLARIZED),
        (ASYMPTOTIC_PER_GENE, SINGLE_ASYMPTOTIC_PER_GENE),
        (ASYMPTOTIC_AGGREGATED, SINGLE_ASYMPTOTIC_AGGREGATED),
        (IMPUTED, SINGLE_IMPUTED),
        (ALPHA_TG, SINGLE_ALPHA_TG),
    ],
    ids=["mk", "polarized", "asymptotic_per_gene", "asymptotic_aggregated", "imputed", "alpha_tg"],
)
def test_single_tsv(result, expected):
    assert format_result(result, TSV) == expected


@pytest.mark.parametrize(
    ("result", "adjusted", "expected"),
    [
        (MK, None, BATCH_MK),
        (MK, ADJUSTED, BATCH_MK_ADJUSTED),
        (POLARIZED, None, BATCH_POLARIZED),
        (POLARIZED, ADJUSTED, BATCH_POLARIZED_ADJUSTED),
        (ASYMPTOTIC_PER_GENE, None, BATCH_ASYMPTOTIC),
        (IMPUTED, None, BATCH_IMPUTED),
        (IMPUTED, ADJUSTED, BATCH_IMPUTED_ADJUSTED),
    ],
    ids=[
        "mk",
        "mk_adjusted",
        "polarized",
        "polarized_adjusted",
        "asymptotic",
        "imputed",
        "imputed_adjusted",
    ],
)
def test_batch_tsv(result, adjusted, expected):
    assert format_batch_results([("g1", result), ("g2", result)], TSV, adjusted) == expected


@pytest.mark.parametrize(
    "result", [MK, ASYMPTOTIC_PER_GENE, IMPUTED], ids=["mk", "asymptotic", "imputed"]
)
def test_single_and_batch_share_one_layout(result):
    """For these types, a batch row is the single row with the gene name in front."""
    single = format_result(result, TSV).split("\n")
    batch = format_batch_results([("g1", result)], TSV).split("\n")
    assert batch == [f"gene\t{single[0]}", f"g1\t{single[1]}"]


@pytest.mark.parametrize(
    ("result", "expected"),
    [(MK, MK.p_value), (POLARIZED, 0.04321), (IMPUTED, 0.0123), (ASYMPTOTIC_PER_GENE, None)],
    ids=["mk", "polarized", "imputed", "asymptotic"],
)
def test_batch_p_value(result, expected):
    assert batch_p_value(result) == expected


def test_batch_p_value_rejects_a_type_with_no_batch_layout():
    with pytest.raises(TypeError, match="AlphaTGResult"):
        batch_p_value(ALPHA_TG)
