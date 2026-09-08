"""Tests for CLI option validation."""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from mkado.cli import app

runner = CliRunner()


class TestOutputPaths:
    @pytest.mark.parametrize("command", ["test", "batch"])
    @pytest.mark.parametrize("option", ["--output", "-O"])
    def test_explicit_stdout_matches_default(self, tmp_path, monkeypatch, command, option):
        monkeypatch.chdir(tmp_path)
        fasta = tmp_path / "test.fa"
        fasta.write_text(">speciesA_1\nATGATGATG\n>speciesA_2\nATGCTGATG\n>speciesB_1\nATGGTGATG\n")
        args = [
            command,
            str(fasta if command == "test" else tmp_path),
            "-i",
            "speciesA",
            "-o",
            "speciesB",
            "--format",
            "json",
        ]
        default = runner.invoke(app, args)
        explicit = runner.invoke(app, [*args, option, "-"])
        assert default.exit_code == 0, default.output
        assert explicit.exit_code == 0, explicit.output
        assert explicit.stdout == default.stdout
        assert not (tmp_path / "-").exists()

    @pytest.mark.parametrize("command", ["test", "batch", "vcf"])
    def test_flag_like_output_is_rejected(self, command):
        result = runner.invoke(app, [command, "--output", "--asymptotic"])
        assert result.exit_code == 2
        assert "looks like a flag" in result.output

    def test_vcf_stdout_reaches_file_validation(self, tmp_path):
        missing = str(tmp_path / "missing.vcf")
        result = runner.invoke(
            app,
            [
                "vcf",
                "--vcf",
                missing,
                "--ref",
                missing,
                "--gff",
                missing,
                "--outgroup-vcf",
                missing,
                "--output",
                "-",
            ],
        )
        assert result.exit_code == 1
        assert "--vcf file not found" in result.output

    def test_other_file_options_still_reject_dash(self):
        result = runner.invoke(app, ["test", "--plot-asymptotic", "-"])
        assert result.exit_code == 2
        assert "looks like a flag" in result.output

    def test_output_file_matches_stdout(self, tmp_path):
        fasta = tmp_path / "test.fa"
        fasta.write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
        args = ["test", str(fasta), "-i", "speciesA", "-o", "speciesB", "-f", "json"]
        default = runner.invoke(app, args)
        target = tmp_path / "results.json"
        saved = runner.invoke(app, [*args, "--output", str(target)])
        assert default.exit_code == saved.exit_code == 0
        assert target.read_text() == default.stdout
        assert saved.stdout == ""


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
        # Should show the singleton exclusion message
        assert "Excluding singletons" in result.output
        assert "min frequency" in result.output

    def test_batch_no_singletons_with_alpha_tg_succeeds(self, tmp_path: Path) -> None:
        """Test that --no-singletons and --alpha-tg can be used together."""
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
                "--no-singletons",
                "--alpha-tg",
            ],
        )

        # Should succeed
        assert result.exit_code == 0


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
