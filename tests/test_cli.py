"""Tests for CLI option validation."""

from pathlib import Path

import pytest
import typer
from typer.testing import CliRunner

from mkado.cli import app, cli_error

runner = CliRunner()


class TestOptionValidation:
    """Tests for CLI option compatibility validation."""

    def test_asymptotic_with_min_freq_error(self, tmp_path: Path) -> None:
        """Test that --asymptotic and --min-freq cannot be used together."""
        # Create a minimal test file
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--min-freq",
                "0.1",
            ],
        )

        assert result.exit_code == 1
        assert "--min-freq cannot be used with --asymptotic" in result.output
        assert "--freq-cutoffs" in result.output

    def test_alpha_tg_with_min_freq_succeeds(self, tmp_path: Path) -> None:
        """Test that --alpha-tg and --min-freq can be used together."""
        # Create a minimal test directory with alignment
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATGATGATGATG
>speciesA_2
ATGCTGATGATGATGATG
>speciesB_1
ATGGTGATGATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--alpha-tg",
                "--min-freq",
                "0.1",
            ],
        )

        # Should succeed
        assert result.exit_code == 0

    def test_alpha_tg_with_asymptotic_error(self, tmp_path: Path) -> None:
        """Test that --alpha-tg and --asymptotic cannot be used together."""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--alpha-tg",
                "--asymptotic",
            ],
        )

        assert result.exit_code == 1
        assert "--alpha-tg and --asymptotic are mutually exclusive" in result.output

    def test_asymptotic_with_polarize_match_error(self, tmp_path: Path) -> None:
        """Test that --asymptotic and --polarize-match cannot be used together."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
>speciesC_1
ATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--polarize-match",
                "speciesC",
            ],
        )

        assert result.exit_code == 1
        assert "Polarized asymptotic test not supported" in result.output

    def test_asymptotic_without_min_freq_succeeds(self, tmp_path: Path) -> None:
        """Test that --asymptotic works without --min-freq."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATGATGATGATG
>speciesA_2
ATGCTGATGATGATGATG
>speciesB_1
ATGGTGATGATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
            ],
        )

        # Should succeed (exit code 0)
        assert result.exit_code == 0

    def test_min_freq_without_asymptotic_succeeds(self, tmp_path: Path) -> None:
        """Test that --min-freq works without --asymptotic (standard MK test)."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesA_2
ATGCTGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--min-freq",
                "0.1",
            ],
        )

        # Should succeed
        assert result.exit_code == 0


class TestBatchOptionValidation:
    """Tests for batch command option validation."""

    def test_batch_asymptotic_with_polarize_match_error(self, tmp_path: Path) -> None:
        """Test that batch --asymptotic and --polarize-match cannot be used together."""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        fasta = alignment_dir / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
>speciesC_1
ATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--polarize-match",
                "speciesC",
            ],
        )

        assert result.exit_code == 1
        assert "--asymptotic and --polarize-match are mutually exclusive" in result.output


class TestNoSingletonsOption:
    """Tests for --no-singletons option."""

    def test_no_singletons_with_min_freq_error(self, tmp_path: Path) -> None:
        """Test that --no-singletons and --min-freq cannot be used together."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--no-singletons",
                "--min-freq",
                "0.1",
            ],
        )

        assert result.exit_code == 1
        assert "--no-singletons and --min-freq cannot be used together" in result.output

    def test_no_singletons_with_asymptotic_error(self, tmp_path: Path) -> None:
        """Test that --no-singletons and --asymptotic cannot be used together."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--no-singletons",
                "--asymptotic",
            ],
        )

        assert result.exit_code == 1
        assert "--no-singletons cannot be used with --asymptotic" in result.output

    def test_no_singletons_succeeds(self, tmp_path: Path) -> None:
        """Test that --no-singletons works correctly."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesA_2
ATGCTGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--no-singletons",
            ],
        )

        assert result.exit_code == 0
        assert "Excluding singletons" in result.output

    def test_batch_no_singletons_drops_the_singleton_for_alpha_tg(self, tmp_path: Path) -> None:
        """The flag reaches the alpha-TG extraction path, which it used not to."""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()
        (alignment_dir / "gene0.fa").write_text(
            ">speciesA_1\nATGGCCAAA\n>speciesA_2\nATGGCCAAA\n"
            ">speciesA_3\nATGGCCAAA\n>speciesA_4\nATGGCTAAA\n"
            ">speciesB_1\nATGGCCAAA\n"
        )
        argv = ["batch", str(alignment_dir), "-i", "speciesA", "-o", "speciesB", "--alpha-tg"]
        kept = runner.invoke(app, argv)
        dropped = runner.invoke(app, [*argv, "--no-singletons"])
        assert kept.exit_code == 0
        assert dropped.exit_code == 0
        assert kept.stdout.splitlines()[1].split("\t")[3] == "1"
        assert dropped.stdout.splitlines()[1].split("\t")[3] == "0"


class TestCodeTable:
    """Tests for --code-table option and mkado codes command."""

    def test_codes_command(self) -> None:
        """Test that mkado codes lists available tables."""
        result = runner.invoke(app, ["codes"])
        assert result.exit_code == 0
        assert "Standard" in result.output
        assert "Vertebrate Mitochondrial" in result.output
        assert "vertebrate-mito" in result.output

    def test_code_table_by_name(self, tmp_path: Path) -> None:
        """Test --code-table accepts a name alias."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATGATGATGATG
>speciesA_2
ATGCTGATGATGATGATG
>speciesB_1
ATGGTGATGATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--code-table",
                "vertebrate-mito",
            ],
        )
        assert result.exit_code == 0

    def test_code_table_by_id(self, tmp_path: Path) -> None:
        """Test --code-table accepts a numeric ID."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATGATGATGATG
>speciesA_2
ATGCTGATGATGATGATG
>speciesB_1
ATGGTGATGATGATGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--code-table",
                "2",
            ],
        )
        assert result.exit_code == 0

    def test_code_table_unknown_name_error(self, tmp_path: Path) -> None:
        """Test --code-table rejects an unknown name."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--code-table",
                "not-a-real-code",
            ],
        )
        assert result.exit_code == 1
        assert "Unknown genetic code" in result.output

    def test_code_table_unknown_id_error(self, tmp_path: Path) -> None:
        """Test --code-table rejects an unknown numeric ID."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(""">speciesA_1
ATGATGATG
>speciesB_1
ATGGTGATG
""")

        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--code-table",
                "99",
            ],
        )
        assert result.exit_code == 1
        assert "Unknown genetic code" in result.output


def _make_asymptotic_alignment_dir(tmp_path: Path) -> Path:
    """Build a tiny multi-gene alignment dir suitable for batch -a."""
    alignment_dir = tmp_path / "alignments"
    alignment_dir.mkdir()
    # Three genes; each has a couple of polymorphisms across frequency
    for i, codons in enumerate(["ATGTTT", "ATGTTC", "ATGTTA"]):
        fa = alignment_dir / f"gene{i}.fa"
        fa.write_text(
            f">speciesA_1\n{codons}\n>speciesA_2\nATG{codons[3:]}\n"
            f">speciesA_3\n{codons}\n>speciesB_1\nATGTTG\n"
        )
    return alignment_dir


class TestCiMethodOption:
    """Tests for --ci-method flag."""

    def _make_alignment_dir(self, tmp_path: Path) -> Path:
        return _make_asymptotic_alignment_dir(tmp_path)

    def test_invalid_ci_method_rejected_in_test(self, tmp_path: Path) -> None:
        fasta = tmp_path / "test.fa"
        fasta.write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
        result = runner.invoke(
            app,
            ["test", str(fasta), "-i", "speciesA", "-o", "speciesB", "--ci-method", "bogus"],
        )
        assert result.exit_code == 1
        assert "Invalid --ci-method" in result.output

    def test_invalid_ci_method_rejected_in_batch(self, tmp_path: Path) -> None:
        alignment_dir = self._make_alignment_dir(tmp_path)
        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--ci-method",
                "bogus",
            ],
        )
        assert result.exit_code == 1
        assert "Invalid --ci-method" in result.output

    def test_ci_method_default_is_monte_carlo_in_tsv(self, tmp_path: Path) -> None:
        """Default --ci-method should produce ci_method=monte-carlo in batch -a output."""
        alignment_dir = self._make_alignment_dir(tmp_path)
        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--format",
                "tsv",
            ],
        )
        assert result.exit_code == 0
        assert "ci_method" in result.output
        assert "monte-carlo" in result.output

    def test_ci_method_bootstrap_appears_in_tsv(self, tmp_path: Path) -> None:
        alignment_dir = self._make_alignment_dir(tmp_path)
        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--ci-method",
                "bootstrap",
                "--bootstrap",
                "20",
                "--format",
                "tsv",
            ],
        )
        assert result.exit_code == 0
        assert "ci_method" in result.output
        assert "bootstrap" in result.output


class TestSfsModeOption:
    """Tests for the --sfs-mode flag (Uricchio et al. 2019 cumulative SFS)."""

    def _make_alignment_dir(self, tmp_path: Path) -> Path:
        return _make_asymptotic_alignment_dir(tmp_path)

    def test_invalid_sfs_mode_rejected_in_test(self, tmp_path: Path) -> None:
        fasta = tmp_path / "test.fa"
        fasta.write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--sfs-mode",
                "bogus",
            ],
        )
        assert result.exit_code == 1
        assert "Invalid --sfs-mode" in result.output

    def test_invalid_sfs_mode_rejected_in_batch(self, tmp_path: Path) -> None:
        alignment_dir = self._make_alignment_dir(tmp_path)
        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--sfs-mode",
                "bogus",
            ],
        )
        assert result.exit_code == 1
        assert "Invalid --sfs-mode" in result.output

    def test_sfs_mode_default_is_at_in_tsv(self, tmp_path: Path) -> None:
        """Default --sfs-mode produces sfs_mode=at in batch -a output."""
        alignment_dir = self._make_alignment_dir(tmp_path)
        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--format",
                "tsv",
            ],
        )
        assert result.exit_code == 0
        assert "sfs_mode" in result.output
        assert "at" in result.output

    def test_sfs_mode_above_appears_in_tsv(self, tmp_path: Path) -> None:
        alignment_dir = self._make_alignment_dir(tmp_path)
        result = runner.invoke(
            app,
            [
                "batch",
                str(alignment_dir),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--sfs-mode",
                "above",
                "--format",
                "tsv",
            ],
        )
        assert result.exit_code == 0
        assert "sfs_mode" in result.output
        assert "above" in result.output

    def test_sfs_mode_above_in_test_command_succeeds(self, tmp_path: Path) -> None:
        fasta = tmp_path / "test.fa"
        fasta.write_text(
            ">speciesA_1\nATGTTTGCAGCAGCAGCAGCAGCA\n"
            ">speciesA_2\nATGTTCGCAGCAGCAGCAGCAGCA\n"
            ">speciesA_3\nATGTTAGCAGCAGCAGCAGCAGCA\n"
            ">speciesA_4\nATGTTTGCAGCAGCAGCAGCAGCA\n"
            ">speciesA_5\nATGTTGGCAGCAGCAGCAGCAGCA\n"
            ">speciesB_1\nATGCTGGCAGCAGCAGCAGCAGCA\n"
        )
        result = runner.invoke(
            app,
            [
                "test",
                str(fasta),
                "-i",
                "speciesA",
                "-o",
                "speciesB",
                "--asymptotic",
                "--sfs-mode",
                "above",
                "--format",
                "tsv",
            ],
        )
        assert result.exit_code == 0
        assert "sfs_mode" in result.output
        assert "above" in result.output


def _run_duplicate_name_batch(tmp_path: Path, output_format: str):
    """Run a batch over two files whose stems collide, so both become gene 'gene1'."""
    alignments = tmp_path / "alignments"
    alignments.mkdir()
    (alignments / "gene1.fa").write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
    (alignments / "gene1.fasta").write_text(">speciesA_1\nATGAAAGGG\n>speciesB_1\nATGCCCGGG\n")
    return runner.invoke(
        app,
        # A single glob can match two extensions, so one directory yields two 'gene1' results.
        [
            "batch",
            str(alignments),
            "--pattern",
            "*.f*",
            "-i",
            "speciesA",
            "-o",
            "speciesB",
            "--format",
            output_format,
        ],
    )


class TestBatchDuplicateGeneNames:
    """Two input files that share a stem give two results with the same gene name."""

    def test_json_reports_the_collision(self, tmp_path: Path) -> None:
        result = _run_duplicate_name_batch(tmp_path, "json")
        assert result.exit_code == 1
        assert "gene1" in result.output
        assert "Traceback" not in result.output

    def test_tsv_keeps_both_rows(self, tmp_path: Path) -> None:
        result = _run_duplicate_name_batch(tmp_path, "tsv")
        assert result.exit_code == 0
        assert [line.split("\t")[0] for line in result.stdout.strip().splitlines()[1:]] == [
            "gene1",
            "gene1",
        ]

    def test_warns_and_names_the_colliding_files(self, tmp_path: Path) -> None:
        """The warning arrives before the run, and names the files the formatter cannot see."""
        result = _run_duplicate_name_batch(tmp_path, "tsv")
        assert "gene1.fa" in result.output
        assert "gene1.fasta" in result.output


class TestCliError:
    def test_prefixes_the_message_and_returns_the_exit(self, capsys: pytest.CaptureFixture) -> None:
        exit_exception = cli_error("the thing went wrong")
        assert isinstance(exit_exception, typer.Exit)
        assert exit_exception.exit_code == 1
        assert capsys.readouterr().err == "Error: the thing went wrong\n"


def test_batch_with_no_results_succeeds(tmp_path: Path) -> None:
    """Unanalysable input warns rather than errors, so an empty result set is not a failure."""
    alignments = tmp_path / "alignments"
    alignments.mkdir()
    (alignments / "gene1.fa").write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
    result = runner.invoke(app, ["batch", str(alignments), "-i", "nomatch", "-o", "alsonomatch"])
    assert result.exit_code == 0
    assert "No results to display" in result.output
