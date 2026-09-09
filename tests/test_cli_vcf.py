"""End-to-end tests for the vcf command."""

from __future__ import annotations

import dataclasses
import json
import logging
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mkado.cli import app

runner = CliRunner()

PER_GENE_HEADER = "gene\tDn\tDs\tPn\tPs\tp_value\tp_value_adjusted\tNI\talpha\tDoS\tLn\tLs\tomega"
SINGLE_HEADER = "Dn\tDs\tPn\tPs\tp_value\tNI\talpha\tDoS\tLn\tLs\tomega"
ASYMPTOTIC_HEADER = "Dn\tDs\tPn\tPs\talpha_asymptotic\tCI_low\tCI_high\tmodel\tnum_genes"
IMPUTED_HEADER = "Dn\tDs\tPn\tPs\tPwd\tPn_neutral\talpha\tp_value\tcutoff"
ALPHA_TG_HEADER = "Dn\tDs\tPn\tPs\talpha_TG\tNI_TG\tCI_low\tCI_high\tnum_genes"
G_PLUS_ROW = "g_plus\t1\t1\t2\t1\t1\t1\t2.000000\t-1.000000\t-0.166667\tNA\tNA\tNA"
G_PLUS_SINGLE_ROW = "1\t1\t2\t1\t1\t2.000000\t-1.000000\t-0.166667\tNA\tNA\tNA"
# Keep the bootstrap small: the defaults multiply into thousands of replicates.
FAST = ["--bootstrap", "5"]

_MKADO_LOGGERS = ("mkado.io.vcf", "mkado.io.plotting")


@pytest.fixture(autouse=True)
def _restore_logging():
    """The vcf command raises logger levels under --verbose; undo that after each test."""
    yield
    for name in _MKADO_LOGGERS:
        logging.getLogger(name).setLevel(logging.NOTSET)


def invoke(dataset, gff: Path, *extra: str, workers: str = "1"):
    args = [
        "vcf",
        "--vcf",
        str(dataset.ingroup_vcf),
        "--ref",
        str(dataset.ref_fasta),
        "--gff",
        str(gff),
        "--outgroup-vcf",
        str(dataset.outgroup_vcf),
        "--workers",
        workers,
        *extra,
    ]
    return runner.invoke(app, args)


def rows(stdout: str) -> dict[str, tuple[int, int, int, int]]:
    """Map gene to (Dn, Ds, Pn, Ps) from per-gene TSV output."""
    lines = stdout.strip().splitlines()
    assert lines[0] == PER_GENE_HEADER
    table = {}
    for line in lines[1:]:
        fields = line.split("\t")
        table[fields[0]] = tuple(int(x) for x in fields[1:5])
    return table


def header_and_row(stdout: str) -> tuple[str, str]:
    lines = stdout.strip().splitlines()
    assert len(lines) == 2, lines
    return lines[0], lines[1]


def expected_rows(genome) -> dict[str, tuple[int, int, int, int]]:
    return {gene_id: expected.counts for gene_id, expected in genome.expected.items()}


def warning_line(genome) -> str:
    return f"Warning: {genome.expected['g_plus'].warning}"


class TestValidation:
    @pytest.mark.parametrize(
        ("extra", "message"),
        [
            (["-f", "xml"], "Error: Invalid format 'xml'."),
            (["--ci-method", "jack"], "Error: Invalid --ci-method 'jack'."),
            (["--sfs-mode", "below"], "Error: Invalid --sfs-mode 'below'."),
            (["-a", "--min-freq", "0.1"], "Error: --min-freq cannot be used with --asymptotic."),
            (["-a", "--no-singletons"], "Error: --no-singletons cannot be used with --asymptotic."),
            (["--imputed", "-a"], "Error: --imputed and --asymptotic are mutually exclusive."),
            (["--alpha-tg", "-a"], "Error: --alpha-tg and --asymptotic are mutually exclusive."),
            (
                ["--imputed", "--alpha-tg"],
                "Error: --imputed and --alpha-tg are mutually exclusive.",
            ),
            (
                ["--imputed", "--no-singletons"],
                "Error: --no-singletons cannot be used with --imputed.",
            ),
            (
                ["--no-singletons", "--min-freq", "0.1"],
                "Error: --no-singletons and --min-freq cannot be used together.",
            ),
            (["--code-table", "99"], "Error: Unknown genetic code table 99."),
            (["--code-table", "martian"], "Error: Unknown genetic code 'martian'."),
            (["--freq-cutoffs", "0.1"], "Error: Invalid frequency cutoffs '0.1'."),
            (["--freq-cutoffs", "a,b"], "Error: Invalid frequency cutoffs 'a,b'."),
        ],
    )
    def test_option_errors(self, genome, gff_chr1, extra, message):
        result = invoke(genome, gff_chr1, *extra)
        assert result.exit_code == 1
        assert message in result.output

    @pytest.mark.parametrize(
        ("option", "field"),
        [
            ("--vcf", "ingroup_vcf"),
            ("--ref", "ref_fasta"),
            ("--gff", None),
            ("--outgroup-vcf", "outgroup_vcf"),
        ],
    )
    def test_missing_file(self, genome, gff_chr1, tmp_path, option, field):
        missing = tmp_path / "missing"
        if field is None:
            result = invoke(genome, missing)
        else:
            result = invoke(dataclasses.replace(genome, **{field: missing}), gff_chr1)
        assert result.exit_code == 1
        assert f"Error: {option} file not found: {missing}" in result.output

    def test_gene_list_looks_like_flag(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--gene-list", "-a")
        assert result.exit_code == 2
        assert "looks like a flag" in result.output


class TestGeneSelection:
    def test_single_gene(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--gene", "g_plus")
        assert result.exit_code == 0
        assert "Found 1 genes in annotation" in result.output
        assert warning_line(genome) in result.output
        assert header_and_row(result.stdout) == (SINGLE_HEADER, G_PLUS_SINGLE_ROW)

    def test_gene_list(self, genome, gff_chr1, tmp_path):
        gene_list = tmp_path / "genes.txt"
        gene_list.write_text("g_plus\n\n# a comment\ng_minus\n")
        result = invoke(genome, gff_chr1, "--gene-list", str(gene_list))
        assert result.exit_code == 0
        assert f"Loaded 2 gene IDs from {gene_list}" in result.output
        assert "Found 2 genes in annotation" in result.output
        assert set(rows(result.stdout)) == {"g_plus", "g_minus"}

    def test_gene_list_missing(self, genome, gff_chr1, tmp_path):
        missing = tmp_path / "nope.txt"
        result = invoke(genome, gff_chr1, "--gene-list", str(missing))
        assert result.exit_code == 1
        assert f"Error: Gene list file not found: {missing}" in result.output

    def test_unknown_gene(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--gene", "nope")
        assert result.exit_code == 1
        assert "Error: No valid CDS regions found in GFF3" in result.output

    def test_no_valid_cds(self, genome, gff_invalid_len):
        result = invoke(genome, gff_invalid_len)
        assert result.exit_code == 1
        assert "Error: No valid CDS regions found in GFF3" in result.output


class TestPerGeneOutput:
    def test_tsv(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1)
        assert result.exit_code == 0
        assert "Found 5 genes in annotation" in result.output
        assert warning_line(genome) in result.output
        assert result.stdout.strip().splitlines()[1] == G_PLUS_ROW
        assert rows(result.stdout) == expected_rows(genome)

    def test_json(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "-f", "json")
        assert result.exit_code == 0
        data = json.loads(result.stdout)
        assert set(data) == set(genome.expected)
        g_plus = data["g_plus"]
        assert (g_plus["dn"], g_plus["ds"], g_plus["pn"], g_plus["ps"]) == (1, 1, 2, 1)
        assert g_plus["ni"] == 2.0
        assert g_plus["p_value_adjusted"] == 1.0

    def test_pretty(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "-f", "pretty")
        assert result.exit_code == 0
        assert "=== g_plus ===" in result.stdout
        assert "  Divergence:    Dn=1, Ds=1" in result.stdout
        assert "  Neutrality Index (NI):  2.0000" in result.stdout
        assert "  p-value (BH adj):     1" in result.stdout

    def test_output_file(self, genome, gff_chr1, tmp_path):
        out = tmp_path / "out.tsv"
        result = invoke(genome, gff_chr1, "--output", str(out))
        assert result.exit_code == 0
        assert f"Results saved to {out}" in result.output
        assert "g_plus" not in result.stdout
        assert out.read_text().splitlines()[1] == G_PLUS_ROW

    def test_output_dash_is_stdout(self, genome, gff_chr1):
        """The --output help text documents '-' as stdout."""
        result = invoke(genome, gff_chr1, "--output", "-")
        assert result.exit_code == 0
        assert "Results saved to" not in result.output
        assert result.stdout.strip().splitlines()[1] == G_PLUS_ROW

    def test_worker_error_reported_and_other_genes_continue(self, genome, gff_bad_chrom):
        result = invoke(genome, gff_bad_chrom)
        assert result.exit_code == 0
        assert "Error processing g_nochrom:" in result.output
        assert "chrZ" in result.output
        assert rows(result.stdout) == {"g_plus": (1, 1, 2, 1)}

    def test_every_gene_failing_exits_nonzero(self, genome, gff_bad_chrom_only):
        """An empty result set is normal, but not when it is empty because every gene failed."""
        result = invoke(genome, gff_bad_chrom_only)
        assert result.exit_code == 1
        assert "Error processing g_nochrom:" in result.output
        assert "No results to display" in result.output

    def test_vcf_with_invalid_content(self, genome, gff_chr1, not_a_vcf):
        result = invoke(dataclasses.replace(genome, ingroup_vcf=not_a_vcf), gff_chr1)
        assert result.exit_code == 1
        assert "Error processing g_plus:" in result.output
        assert "No results to display" in result.output


class TestSingleGeneModes:
    def test_json(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--gene", "g_plus", "-f", "json")
        assert result.exit_code == 0
        assert json.loads(result.stdout)["ni"] == 2.0

    def test_pretty(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--gene", "g_plus", "-f", "pretty")
        assert result.exit_code == 0
        assert "  Neutrality Index (NI):  2.0000" in result.stdout

    @pytest.mark.parametrize(
        ("extra", "header"),
        [
            (["-a"], ASYMPTOTIC_HEADER),
            (["--imputed"], IMPUTED_HEADER),
            (["--alpha-tg"], ALPHA_TG_HEADER),
        ],
    )
    def test_single_gene_runs_requested_mode(self, genome, gff_chr1, extra, header):
        """A mode asked for alongside --gene runs that mode on the one selected gene."""
        result = invoke(genome, gff_chr1, "--gene", "g_plus", *extra, *FAST)
        assert result.exit_code == 0
        # Only the aggregate branches announce their gene count, so this pins the path.
        assert "Using 1 genes for" in result.output
        got_header, row = header_and_row(result.stdout)
        assert got_header.startswith(header)
        assert row.split("\t")[:4] == ["1", "1", "2", "1"]

    def test_asymptotic_per_gene(self, genome, gff_chr1):
        result = invoke(
            genome, gff_chr1, "--gene", "g_plus", "-a", "--per-gene", "--bootstrap", "2"
        )
        assert result.exit_code == 0
        header, row = header_and_row(result.stdout)
        assert header.startswith(ASYMPTOTIC_HEADER)
        assert row.startswith("1\t1\t2\t1\t-1.000000\t-1.000000\t-1.000000\texponential\t1\t")
        assert row.endswith("\tmonte-carlo\tat")

    def test_imputed_per_gene(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--gene", "g_plus", "--imputed", "--per-gene", *FAST)
        assert result.exit_code == 0
        header, row = header_and_row(result.stdout)
        assert header.startswith(IMPUTED_HEADER)
        assert row.startswith("1\t1\t2\t1\t0.00\t2.00\t-1.000000\t1\t0.15\t")


class TestAggregateModes:
    def test_asymptotic(self, genome, gff_all):
        result = invoke(genome, gff_all, "-a", *FAST)
        assert result.exit_code == 0
        assert "Found 17 genes in annotation" in result.output
        assert "Using 17 genes for aggregated asymptotic" in result.output
        header, row = header_and_row(result.stdout)
        assert header.startswith(ASYMPTOTIC_HEADER)
        fields = row.split("\t")
        assert fields[:4] == ["15", "3", "10", "16"]
        assert fields[8] == "17"
        assert fields[-2:] == ["monte-carlo", "at"]

    def test_asymptotic_ci_and_sfs_options(self, genome, gff_all):
        result = invoke(
            genome, gff_all, "-a", *FAST, "--ci-method", "bootstrap", "--sfs-mode", "above"
        )
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        assert row.endswith("\tbootstrap\tabove")

    def test_asymptotic_freq_cutoffs(self, genome, gff_all):
        result = invoke(genome, gff_all, "-a", *FAST, "--freq-cutoffs", "0.2,0.8")
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        assert row.split("\t")[:4] == ["15", "3", "10", "16"]

    def test_imputed(self, genome, gff_all):
        result = invoke(genome, gff_all, "--imputed", *FAST)
        assert result.exit_code == 0
        assert "Using 17 genes for aggregated imputed" in result.output
        header, row = header_and_row(result.stdout)
        assert header.startswith(IMPUTED_HEADER)
        assert row.startswith("15\t3\t10\t16\t0.00\t10.00\t0.875000\t")
        assert row.split("\t")[8] == "0.15"

    def test_imputed_min_freq_sets_cutoff(self, genome, gff_all):
        result = invoke(genome, gff_all, "--imputed", *FAST, "--min-freq", "0.3")
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        fields = row.split("\t")
        assert fields[:4] == ["15", "3", "2", "1"]
        assert fields[8] == "0.3"

    def test_alpha_tg(self, genome, gff_all):
        result = invoke(genome, gff_all, "--alpha-tg", *FAST)
        assert result.exit_code == 0
        assert "Using 17 genes for alpha-TG" in result.output
        header, row = header_and_row(result.stdout)
        assert header.startswith(ALPHA_TG_HEADER)
        assert row.startswith("15\t3\t10\t16\t0.821429\t0.178571\t")
        assert row.split("\t")[8] == "17"

    def test_asymptotic_per_gene_table(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "-a", "--per-gene", "--bootstrap", "2")
        assert result.exit_code == 0
        lines = result.stdout.strip().splitlines()
        assert lines[0].startswith("gene\tDn\tDs\talpha_asymptotic\tCI_low\tCI_high\tmodel")
        assert lines[1].startswith("g_plus\t1\t1\t-1.000000\t-1.000000\t-1.000000\texponential\t")
        assert len(lines) == 6

    def test_asymptotic_per_gene_honors_freq_cutoffs(self, genome, gff_chr1, fitter_calls):
        result = invoke(
            genome, gff_chr1, "-a", "--per-gene", "--bootstrap", "2", "--freq-cutoffs", "0.2,0.8"
        )
        assert result.exit_code == 0
        assert [call["frequency_cutoffs"] for call in fitter_calls] == [(0.2, 0.8)] * 5

    def test_imputed_per_gene_table(self, genome, gff_chr1):
        """A per-gene batch carries a real adjusted p-value, spliced after p_value."""
        result = invoke(genome, gff_chr1, "--imputed", "--per-gene", *FAST)
        assert result.exit_code == 0
        lines = result.stdout.strip().splitlines()
        assert lines[0].startswith(
            "gene\t"
            + IMPUTED_HEADER.replace("p_value\tcutoff", "p_value\tp_value_adjusted\tcutoff")
        )
        assert len(lines) == 6
        fields = lines[1].split("\t")
        assert fields[:9] == ["g_plus", "1", "1", "2", "1", "0.00", "2.00", "-1.000000", "1"]
        float(fields[9])  # p_value_adjusted is a real number, not the old hardcoded 1.0
        assert fields[10] == "0.15"


class TestPlots:
    def test_volcano_saved(self, genome, gff_chr1, tmp_path):
        plot = tmp_path / "volcano.png"
        result = invoke(genome, gff_chr1, "--volcano", str(plot))
        assert result.exit_code == 0
        assert f"Volcano plot saved to {plot}" in result.output
        assert plot.stat().st_size > 0

    def test_volcano_without_valid_points(self, genome, gff_tiny, tmp_path):
        plot = tmp_path / "volcano.png"
        result = invoke(genome, gff_tiny, "--volcano", str(plot))
        assert result.exit_code == 0
        assert "Could not generate volcano plot:" in result.output
        assert not plot.exists()

    def test_asymptotic_plot_without_bins(self, genome, gff_tiny, tmp_path):
        plot = tmp_path / "asymptotic.png"
        result = invoke(genome, gff_tiny, "-a", *FAST, "--plot-asymptotic", str(plot))
        assert result.exit_code == 0
        assert "Could not generate asymptotic plot:" in result.output
        assert not plot.exists()

    def test_asymptotic_plot_requires_flag(self, genome, gff_chr1, tmp_path):
        plot = tmp_path / "asymptotic.png"
        result = invoke(genome, gff_chr1, "--plot-asymptotic", str(plot))
        assert result.exit_code == 0
        assert "Warning: --plot-asymptotic requires --asymptotic/-a flag (ignored)" in result.output
        assert not plot.exists()


class TestOptions:
    @pytest.mark.parametrize("table", ["2", "vertebrate-mito"])
    def test_code_table(self, genome, gff_chr1, table):
        result = invoke(genome, gff_chr1, "--code-table", table)
        assert result.exit_code == 0
        table_rows = rows(result.stdout)
        assert table_rows["g_minus"] == (0, 2, 1, 1)
        assert table_rows["g_plus"] == (1, 1, 2, 1)

    def test_min_freq(self, genome, gff_chr1):
        result = invoke(genome, gff_chr1, "--min-freq", "0.2")
        assert result.exit_code == 0
        assert rows(result.stdout)["g_plus"] == (1, 1, 2, 0)

    def test_verbose_enables_debug_logging(self, genome, gff_chr1, ingroup_vcf_plain, caplog):
        """An unindexed VCF makes every region query fail; --verbose surfaces the reason."""
        caplog.handler.setLevel(logging.DEBUG)
        dataset = dataclasses.replace(genome, ingroup_vcf=ingroup_vcf_plain)
        result = invoke(dataset, gff_chr1, "--gene", "g_plus", "--verbose")
        assert result.exit_code == 0
        messages = [r.message for r in caplog.records if r.name == "mkado.io.vcf"]
        assert any("ingroup VCF query failed for g_plus" in m for m in messages)

    def test_without_verbose_no_debug_logging(self, genome, gff_chr1, ingroup_vcf_plain, caplog):
        caplog.handler.setLevel(logging.DEBUG)
        dataset = dataclasses.replace(genome, ingroup_vcf=ingroup_vcf_plain)
        result = invoke(dataset, gff_chr1, "--gene", "g_plus")
        assert result.exit_code == 0
        assert not [r for r in caplog.records if r.name == "mkado.io.vcf"]

    def test_bgzipped_reference(self, genome, gff_chr1, bgzipped_ref):
        result = invoke(dataclasses.replace(genome, ref_fasta=bgzipped_ref), gff_chr1)
        assert result.exit_code == 0
        assert rows(result.stdout) == expected_rows(genome)


class TestParallel:
    def test_two_workers_match_sequential(self, genome, gff_all):
        sequential = invoke(genome, gff_all)
        parallel = invoke(genome, gff_all, workers="2")
        assert parallel.exit_code == 0
        assert warning_line(genome) in parallel.output
        assert rows(parallel.stdout) == rows(sequential.stdout)
        assert len(rows(parallel.stdout)) == 17

    def test_gene_error_inside_chunk(self, genome, gff_all_with_bad_chrom):
        result = invoke(genome, gff_all_with_bad_chrom, workers="2")
        assert result.exit_code == 0
        assert "Error processing g_nochrom:" in result.output
        assert len(rows(result.stdout)) == 17

    def test_chunk_failure_reports_every_gene(self, genome, gff_all, not_a_vcf):
        dataset = dataclasses.replace(genome, ingroup_vcf=not_a_vcf)
        result = invoke(dataset, gff_all, workers="2")
        assert result.exit_code == 1
        assert "Error processing g_plus:" in result.output
        assert "Error processing tiny11:" in result.output
        assert "No results to display" in result.output


@pytest.fixture(scope="module")
def example_per_gene(example_dataset):
    """One sequential per-gene run on the example data, shared by the tests that need it."""
    return invoke(example_dataset, example_dataset.gff)


class TestExampleData:
    def test_per_gene(self, example_per_gene):
        assert example_per_gene.exit_code == 0
        assert "Found 42 genes in annotation" in example_per_gene.output
        table = rows(example_per_gene.stdout)
        assert len(table) == 42
        assert table["LOC_00000001"] == (0, 0, 5, 2)

    def test_single_gene(self, example_dataset):
        result = invoke(example_dataset, example_dataset.gff, "--gene", "LOC_00000001")
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        assert row == "0\t0\t5\t2\t1\tNA\tNA\t-0.714286\tNA\tNA\tNA"

    def test_asymptotic(self, example_dataset):
        result = invoke(example_dataset, example_dataset.gff, "-a", "--bootstrap", "10")
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        assert row.startswith("20\t36\t9\t11\t1.000000\t1.000000\t1.000000\texponential\t42\t")

    def test_asymptotic_plot_saved(self, example_dataset, tmp_path):
        plot = tmp_path / "asymptotic.png"
        result = invoke(
            example_dataset, example_dataset.gff, "-a", *FAST, "--plot-asymptotic", str(plot)
        )
        assert result.exit_code == 0
        assert f"Asymptotic plot saved to {plot}" in result.output
        assert plot.stat().st_size > 0

    def test_imputed(self, example_dataset):
        result = invoke(example_dataset, example_dataset.gff, "--imputed", *FAST)
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        assert row.startswith("20\t36\t9\t11\t1.67\t7.33\t-0.200000\t1\t0.15\t")

    def test_alpha_tg(self, example_dataset):
        result = invoke(example_dataset, example_dataset.gff, "--alpha-tg", *FAST)
        assert result.exit_code == 0
        _, row = header_and_row(result.stdout)
        assert row.startswith("20\t36\t9\t11\t0.710345\t0.289655\t")

    def test_parallel_matches_sequential(self, example_dataset, example_per_gene):
        parallel = invoke(example_dataset, example_dataset.gff, workers="2")
        assert parallel.exit_code == 0
        assert rows(parallel.stdout) == rows(example_per_gene.stdout)
        assert len(rows(parallel.stdout)) == 42
