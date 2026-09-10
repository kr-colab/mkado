"""Tests for the modes of the batch command."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import Result
from typer.testing import CliRunner

from mkado.cli import app, find_partner_file
from tests.builders import TWO_GENES, TwoGeneDirs, write_two_genes

runner = CliRunner()

# Fewest replicates and bins that still let every fit run.
FAST_ASYMPTOTIC = ["--bootstrap", "5", "--bins", "3"]
LAYOUTS = ["combined", "separate"]
AGGREGATE_MODES = [
    (["--imputed", "--bootstrap", "5"], "aggregated imputed"),
    (["--asymptotic", *FAST_ASYMPTOTIC], "aggregated asymptotic"),
    (["--alpha-tg", "--bootstrap", "5"], "alpha-TG"),
]

# The two genes with every polymorphism removed, so an aggregated asymptotic
# run has no frequency bin to plot.
NO_POLYMORPHISM = {
    gene: {"ingroup": [groups["ingroup"][0]] * 2, "outgroup": groups["outgroup"]}
    for gene, groups in TWO_GENES.items()
}
# The two genes with outgroups one codon short, so every gene fails to align.
SHORT_OUTGROUPS = {
    gene: {**groups, "outgroup": [seq[:-3] for seq in groups["outgroup"]]}
    for gene, groups in TWO_GENES.items()
}


def batch(*args: str | Path) -> Result:
    return runner.invoke(app, ["batch", *map(str, args), "-w", "1"])


def batch_in(dirs: TwoGeneDirs, layout: str, *args: str | Path) -> Result:
    """Run batch on one layout of a two-gene set."""
    if layout == "combined":
        return batch(dirs.combined, "-i", "speciesA", "-o", "speciesB", *args)
    return batch(dirs.separate, *args)


def rows(stdout: str) -> list[list[str]]:
    """The TSV table on stdout, split into cells."""
    return [line.split("\t") for line in stdout.strip().splitlines()]


def relocated(dirs: TwoGeneDirs, root: Path, ingroup_name: str, outgroup_name: str) -> Path:
    """Move gene g1's pair into a fresh directory under new file names."""
    layout = root / "relocated"
    layout.mkdir()
    (dirs.separate / "g1_ingroup.fa").rename(layout / ingroup_name)
    (dirs.separate / "g1_outgroup.fa").rename(layout / outgroup_name)
    return layout


class TestSeparateFiles:
    """The mode used when -i is absent: one ingroup file and one outgroup file per gene."""

    def test_standard_counts(self, two_genes: TwoGeneDirs) -> None:
        result = batch(two_genes.separate)
        assert result.exit_code == 0
        assert "Found 2 ingroup files" in result.output
        table = rows(result.stdout)
        assert table[0][:5] == ["gene", "Dn", "Ds", "Pn", "Ps"]
        assert table[1][:5] == ["g1_ingroup", "1", "1", "0", "2"]
        assert table[2][:5] == ["g2_ingroup", "0", "1", "1", "1"]

    def test_polarize_pattern(self, two_genes: TwoGeneDirs) -> None:
        """The second outgroup matches the ingroup, so every difference is on the outgroup lineage."""
        result = batch(two_genes.separate, "--polarize-pattern", "*_outgroup2.fa")
        assert result.exit_code == 0
        table = rows(result.stdout)
        assert table[0][:7] == [
            "gene",
            "Dn_ingroup",
            "Ds_ingroup",
            "Pn_ingroup",
            "Ps_ingroup",
            "Dn_outgroup",
            "Ds_outgroup",
        ]
        assert table[1][:7] == ["g1_ingroup", "0", "0", "0", "2", "1", "1"]
        assert table[2][:7] == ["g2_ingroup", "0", "0", "1", "1", "0", "1"]

    def test_in_and_out_markers(self, two_genes: TwoGeneDirs, tmp_path: Path) -> None:
        """An _in file pairs with the _out file of the same name."""
        layout = relocated(two_genes, tmp_path, "g1_in.fa", "g1_out.fa")
        result = batch(layout, "--ingroup-pattern", "*_in.fa", "--outgroup-pattern", "*_out.fa")
        assert result.exit_code == 0
        assert rows(result.stdout)[1][:5] == ["g1_in", "1", "1", "0", "2"]

    def test_unmarked_ingroup_pairs_with_stem_outgroup(
        self, two_genes: TwoGeneDirs, tmp_path: Path
    ) -> None:
        """An ingroup file without a marker pairs with <stem>_outgroup.fa."""
        layout = relocated(two_genes, tmp_path, "g1.fa", "g1_outgroup.fa")
        result = batch(layout, "--ingroup-pattern", "g?.fa")
        assert result.exit_code == 0
        assert rows(result.stdout)[1][:5] == ["g1", "1", "1", "0", "2"]

    def test_outgroup_found_by_pattern_when_the_conventional_name_is_absent(
        self, two_genes: TwoGeneDirs, tmp_path: Path
    ) -> None:
        layout = relocated(two_genes, tmp_path, "g1_ingroup.fa", "g1_out.fa")
        result = batch(layout, "--outgroup-pattern", "*_out.fa")
        assert result.exit_code == 0
        assert rows(result.stdout)[1][:5] == ["g1_ingroup", "1", "1", "0", "2"]

    def test_ingroup_without_outgroup_is_skipped_with_a_warning(
        self, two_genes: TwoGeneDirs
    ) -> None:
        (two_genes.separate / "g2_outgroup.fa").unlink()
        result = batch(two_genes.separate)
        assert result.exit_code == 0
        assert "Warning: No outgroup for g2_ingroup.fa" in result.output
        assert [row[0] for row in rows(result.stdout)[1:]] == ["g1_ingroup"]

    def test_missing_second_outgroup_is_skipped_with_a_warning(
        self, two_genes: TwoGeneDirs
    ) -> None:
        (two_genes.separate / "g2_outgroup2.fa").unlink()
        result = batch(two_genes.separate, "--polarize-pattern", "*_outgroup2.fa")
        assert result.exit_code == 0
        assert "Warning: No second outgroup for g2_ingroup.fa" in result.output
        assert [row[0] for row in rows(result.stdout)[1:]] == ["g1_ingroup"]

    def test_no_ingroup_files_is_an_error(self, tmp_path: Path) -> None:
        result = batch(tmp_path)
        assert result.exit_code == 1
        assert "No files matching '*_ingroup.fa'" in result.output

    def test_no_valid_pairs_is_an_error(self, two_genes: TwoGeneDirs) -> None:
        (two_genes.separate / "g1_outgroup.fa").unlink()
        (two_genes.separate / "g2_outgroup.fa").unlink()
        result = batch(two_genes.separate)
        assert result.exit_code == 1
        assert "No valid file pairs found" in result.output

    @pytest.mark.parametrize(
        "flags, first_row",
        [
            (["--imputed", "--per-gene"], ["g1_ingroup", "1", "1", "0", "2"]),
            (["--asymptotic", "--per-gene", *FAST_ASYMPTOTIC], ["g1_ingroup", "1", "1"]),
        ],
    )
    def test_per_gene_modes_give_one_row_per_gene(
        self, two_genes: TwoGeneDirs, flags: list[str], first_row: list[str]
    ) -> None:
        result = batch(two_genes.separate, *flags)
        assert result.exit_code == 0
        table = rows(result.stdout)
        assert len(table) == 3
        assert table[1][: len(first_row)] == first_row


class TestCombinedFiles:
    """The mode used with -i and -o: every species in one file per gene."""

    def test_ingroup_match_requires_outgroup_match(self, two_genes: TwoGeneDirs) -> None:
        result = batch(two_genes.combined, "-i", "speciesA")
        assert result.exit_code == 1
        assert "-o/--outgroup-match required with -i/--ingroup-match" in result.output

    def test_pattern_matching_nothing_is_an_error(self, two_genes: TwoGeneDirs) -> None:
        result = batch(
            two_genes.combined, "-i", "speciesA", "-o", "speciesB", "--pattern", "*.nothere"
        )
        assert result.exit_code == 1
        assert "No files found in" in result.output
        assert "*.nothere" in result.output


@pytest.mark.parametrize("layout", LAYOUTS)
class TestAggregateModes:
    """Modes that pool every gene into one result, in both layouts."""

    @pytest.mark.parametrize("flags, label", AGGREGATE_MODES)
    def test_pools_both_genes(
        self, two_genes: TwoGeneDirs, layout: str, flags: list[str], label: str
    ) -> None:
        result = batch_in(two_genes, layout, *flags)
        assert result.exit_code == 0
        assert f"Using 2 genes for {label}" in result.output
        table = rows(result.stdout)
        assert len(table) == 2
        assert table[1][:4] == ["1", "2", "1", "3"]

    @pytest.mark.parametrize("flags, label", AGGREGATE_MODES)
    def test_every_gene_failing_leaves_no_data(
        self, tmp_path: Path, layout: str, flags: list[str], label: str
    ) -> None:
        """Each failure is reported, and the run ends with no result but no error."""
        dirs = write_two_genes(tmp_path, SHORT_OUTGROUPS)
        result = batch_in(dirs, layout, *flags)
        assert result.exit_code == 0
        assert "Error processing g1" in result.output
        assert "Error processing g2" in result.output
        assert f"Using 0 genes for {label}" in result.output
        assert "No valid gene data extracted" in result.output
        assert result.stdout.strip() == ""

    def test_plot_asymptotic_is_saved(
        self, two_genes: TwoGeneDirs, layout: str, tmp_path: Path
    ) -> None:
        plot = tmp_path / "alpha.png"
        result = batch_in(
            two_genes, layout, "--asymptotic", *FAST_ASYMPTOTIC, "--plot-asymptotic", plot
        )
        assert result.exit_code == 0
        assert f"Asymptotic plot saved to {plot}" in result.output
        assert plot.stat().st_size > 0

    def test_plot_asymptotic_without_bins_reports_and_continues(
        self, tmp_path: Path, layout: str
    ) -> None:
        """With no polymorphism there is no frequency bin to plot, but the run still succeeds."""
        dirs = write_two_genes(tmp_path, NO_POLYMORPHISM)
        plot = tmp_path / "alpha.png"
        result = batch_in(dirs, layout, "--asymptotic", *FAST_ASYMPTOTIC, "--plot-asymptotic", plot)
        assert result.exit_code == 0
        assert "Could not generate plot: No frequency bin data available" in result.output
        assert not plot.exists()


@pytest.mark.parametrize("layout", LAYOUTS)
def test_per_gene_mode_with_every_gene_failing_is_an_error(tmp_path: Path, layout: str) -> None:
    """Per-gene mode reports each failure and exits nonzero when nothing succeeded."""
    dirs = write_two_genes(tmp_path, SHORT_OUTGROUPS)
    result = batch_in(dirs, layout)
    assert result.exit_code == 1
    assert "Error processing g1" in result.output
    assert "No results to display; every gene failed" in result.output


class TestFindPartnerFile:
    """Which file in a separate-files directory belongs to an ingroup file."""

    def test_prefix_does_not_claim_another_genes_file(self, tmp_path: Path) -> None:
        """g1 must not take g10's file, whichever one the directory lists first."""
        candidates = [tmp_path / "g10_outgroup2.fa", tmp_path / "g1_outgroup2.fa"]
        assert find_partner_file(tmp_path / "g1_ingroup.fa", candidates) == (
            tmp_path / "g1_outgroup2.fa"
        )
        assert find_partner_file(tmp_path / "g10_ingroup.fa", candidates) == (
            tmp_path / "g10_outgroup2.fa"
        )

    def test_gene_name_keeps_every_token_before_the_marker(self, tmp_path: Path) -> None:
        """Adh_dmel must not take Adh_dsim's file."""
        candidates = [tmp_path / "Adh_dsim_outgroup2.fa", tmp_path / "Adh_dmel_outgroup2.fa"]
        assert find_partner_file(tmp_path / "Adh_dmel_ingroup.fa", candidates) == (
            tmp_path / "Adh_dmel_outgroup2.fa"
        )

    @pytest.mark.parametrize("ingroup", ["g1_ingroup.fa", "g1_in.fa", "g1.fa"])
    def test_marker_is_stripped_from_the_gene_name(self, tmp_path: Path, ingroup: str) -> None:
        candidates = [tmp_path / "g1_out.fa"]
        assert find_partner_file(tmp_path / ingroup, candidates) == tmp_path / "g1_out.fa"

    def test_no_candidate_for_the_gene(self, tmp_path: Path) -> None:
        assert find_partner_file(tmp_path / "g1_ingroup.fa", [tmp_path / "g2_out.fa"]) is None
