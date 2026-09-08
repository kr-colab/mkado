"""Tests for parallel batch processing."""

from pathlib import Path
from unittest.mock import patch

from mkado.batch_workers import BatchTask, WorkerResult, process_gene
from mkado.cli import get_worker_count, run_parallel_batch


class TestGetWorkerCount:
    """Tests for the get_worker_count helper function."""

    def test_sequential_when_requested(self) -> None:
        """Test that workers=1 always returns 1."""
        assert get_worker_count(1, 100) == 1
        assert get_worker_count(1, 1000) == 1

    def test_sequential_for_few_tasks(self) -> None:
        """Test that small task counts return 1 worker."""
        assert get_worker_count(0, 5) == 1
        assert get_worker_count(0, 9) == 1
        assert get_worker_count(4, 5) == 1

    def test_auto_detection(self) -> None:
        """Test auto-detection of worker count."""
        with patch("os.cpu_count", return_value=8):
            # Auto with many tasks should use cpu_count - 1
            result = get_worker_count(0, 100)
            assert result == 7

    def test_auto_capped_by_task_count(self) -> None:
        """Test that auto mode doesn't exceed task count."""
        with patch("os.cpu_count", return_value=16):
            # With 12 tasks, should use min(15, 12) = 12
            result = get_worker_count(0, 12)
            assert result == 12

    def test_requested_count_capped(self) -> None:
        """Test that requested workers are capped at CPU count."""
        with patch("os.cpu_count", return_value=4):
            assert get_worker_count(8, 100) == 4
            assert get_worker_count(4, 100) == 4
            assert get_worker_count(2, 100) == 2

    def test_handles_no_cpu_count(self) -> None:
        """Test handling when os.cpu_count returns None."""
        with patch("os.cpu_count", return_value=None):
            # Should default to 4 CPUs
            result = get_worker_count(0, 100)
            assert result == 3  # cpu_count - 1 = 4 - 1 = 3


class TestProcessGene:
    """Tests for the process_gene worker function."""

    def test_combined_mode_standard_mk(self, tmp_path: Path) -> None:
        """Test process_gene with combined file mode and standard MK test."""
        # Create test alignment file
        alignment = tmp_path / "gene1.fa"
        alignment.write_text(""">gene1_speciesA_1
ATGATGATG
>gene1_speciesA_2
ATGCTGATG
>gene1_speciesB_1
ATGGTGATG
""")

        task = BatchTask(
            file_path=alignment,
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            reading_frame=1,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1"
        assert result.error is None
        assert result.warning is None
        assert result.result is not None

    def test_combined_mode_no_ingroup(self, tmp_path: Path) -> None:
        """Test process_gene returns warning when no ingroup matches."""
        alignment = tmp_path / "gene1.fa"
        alignment.write_text(""">gene1_speciesB_1
ATGATGATG
""")

        task = BatchTask(
            file_path=alignment,
            ingroup_match="speciesA",  # No matches
            outgroup_match="speciesB",
            reading_frame=1,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1"
        assert result.warning is not None
        assert "No ingroup" in result.warning

    def test_combined_mode_no_outgroup(self, tmp_path: Path) -> None:
        """Test process_gene returns warning when no outgroup matches."""
        alignment = tmp_path / "gene1.fa"
        alignment.write_text(""">gene1_speciesA_1
ATGATGATG
""")

        task = BatchTask(
            file_path=alignment,
            ingroup_match="speciesA",
            outgroup_match="speciesC",  # No matches
            reading_frame=1,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1"
        assert result.warning is not None
        assert "No outgroup" in result.warning

    def test_separate_files_mode(self, tmp_path: Path) -> None:
        """Test process_gene with separate files mode."""
        ingroup = tmp_path / "gene1_ingroup.fa"
        ingroup.write_text(""">seq1
ATGATGATG
>seq2
ATGCTGATG
""")

        outgroup = tmp_path / "gene1_outgroup.fa"
        outgroup.write_text(""">out1
ATGGTGATG
""")

        task = BatchTask(
            file_path=ingroup,
            outgroup_file=outgroup,
            reading_frame=1,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1_ingroup"
        assert result.error is None
        assert result.warning is None
        assert result.result is not None

    def test_separate_files_no_outgroup(self, tmp_path: Path) -> None:
        """Test process_gene returns warning when outgroup file is missing."""
        ingroup = tmp_path / "gene1_ingroup.fa"
        ingroup.write_text(""">seq1
ATGATGATG
""")

        task = BatchTask(
            file_path=ingroup,
            outgroup_file=None,  # Missing outgroup
            reading_frame=1,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1_ingroup"
        assert result.warning is not None
        assert "No outgroup" in result.warning

    def test_extract_only_mode(self, tmp_path: Path) -> None:
        """Test process_gene with extract_only for aggregated asymptotic."""
        alignment = tmp_path / "gene1.fa"
        alignment.write_text(""">gene1_speciesA_1
ATGATGATG
>gene1_speciesA_2
ATGCTGATG
>gene1_speciesB_1
ATGGTGATG
""")

        task = BatchTask(
            file_path=alignment,
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            reading_frame=1,
            extract_only=True,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1"
        assert result.error is None
        assert result.result is not None
        # Result should be PolymorphismData
        assert hasattr(result.result, "polymorphisms")
        assert hasattr(result.result, "dn")
        assert hasattr(result.result, "ds")

    def test_asymptotic_mode(self, tmp_path: Path) -> None:
        """Test process_gene with asymptotic MK test."""
        # Need more sequences for meaningful asymptotic test
        alignment = tmp_path / "gene1.fa"
        alignment.write_text(""">gene1_speciesA_1
ATGATGATGATGATGATG
>gene1_speciesA_2
ATGCTGATGATGATGATG
>gene1_speciesA_3
ATGATGATGATGATGATG
>gene1_speciesA_4
ATGATGATGATGATGATG
>gene1_speciesB_1
ATGGTGATGATGATGATG
>gene1_speciesB_2
ATGGTGATGATGATGATG
""")

        task = BatchTask(
            file_path=alignment,
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            reading_frame=1,
            use_asymptotic=True,
            bins=5,
            bootstrap=10,
        )

        result = process_gene(task)

        assert result.gene_id == "gene1"
        assert result.error is None
        assert result.result is not None
        # Result should be AsymptoticMKResult
        assert hasattr(result.result, "alpha_asymptotic")


class TestRunParallelBatch:
    """Tests for the run_parallel_batch function."""

    def test_sequential_mode(self, tmp_path: Path) -> None:
        """Test batch processing in sequential mode (workers=1)."""
        # Create test files
        for i in range(3):
            f = tmp_path / f"gene{i}.fa"
            f.write_text(f""">gene{i}_speciesA_1
ATGATGATG
>gene{i}_speciesA_2
ATGCTGATG
>gene{i}_speciesB_1
ATGGTGATG
""")

        tasks = [
            BatchTask(
                file_path=tmp_path / f"gene{i}.fa",
                ingroup_match="speciesA",
                outgroup_match="speciesB",
                reading_frame=1,
            )
            for i in range(3)
        ]

        results, warnings, _ = run_parallel_batch(tasks, 1, "Testing")

        assert len(results) == 3
        assert len(warnings) == 0
        gene_ids = {r.gene_id for r in results}
        assert gene_ids == {"gene0", "gene1", "gene2"}

    def test_parallel_mode(self, tmp_path: Path) -> None:
        """Test batch processing in parallel mode."""
        # Create test files
        for i in range(5):
            f = tmp_path / f"gene{i}.fa"
            f.write_text(f""">gene{i}_speciesA_1
ATGATGATG
>gene{i}_speciesA_2
ATGCTGATG
>gene{i}_speciesB_1
ATGGTGATG
""")

        tasks = [
            BatchTask(
                file_path=tmp_path / f"gene{i}.fa",
                ingroup_match="speciesA",
                outgroup_match="speciesB",
                reading_frame=1,
            )
            for i in range(5)
        ]

        results, warnings, _ = run_parallel_batch(tasks, 2, "Testing parallel")

        assert len(results) == 5
        assert len(warnings) == 0
        gene_ids = {r.gene_id for r in results}
        assert gene_ids == {"gene0", "gene1", "gene2", "gene3", "gene4"}

    def test_handles_warnings(self, tmp_path: Path) -> None:
        """Test that warnings are collected properly."""
        # Create one valid file and one that will produce a warning
        valid = tmp_path / "valid.fa"
        valid.write_text(""">valid_speciesA_1
ATGATGATG
>valid_speciesB_1
ATGGTGATG
""")

        # File with only outgroup (will warn about no ingroup)
        no_ingroup = tmp_path / "no_ingroup.fa"
        no_ingroup.write_text(""">no_ingroup_speciesB_1
ATGATGATG
""")

        tasks = [
            BatchTask(
                file_path=valid,
                ingroup_match="speciesA",
                outgroup_match="speciesB",
                reading_frame=1,
            ),
            BatchTask(
                file_path=no_ingroup,
                ingroup_match="speciesA",
                outgroup_match="speciesB",
                reading_frame=1,
            ),
        ]

        results, warnings, _ = run_parallel_batch(tasks, 1, "Testing warnings")

        assert len(results) == 1
        assert len(warnings) == 1
        assert "No ingroup" in warnings[0]

    def test_same_results_sequential_vs_parallel(self, tmp_path: Path) -> None:
        """Test that sequential and parallel modes produce the same results."""
        # Create test files with deterministic content
        for i in range(4):
            f = tmp_path / f"gene{i}.fa"
            f.write_text(f""">gene{i}_speciesA_1
ATGATGATGATG
>gene{i}_speciesA_2
ATGCTGATGATG
>gene{i}_speciesB_1
ATGGTGATGATG
""")

        tasks = [
            BatchTask(
                file_path=tmp_path / f"gene{i}.fa",
                ingroup_match="speciesA",
                outgroup_match="speciesB",
                reading_frame=1,
            )
            for i in range(4)
        ]

        # Run in sequential mode
        seq_results, seq_warnings, _ = run_parallel_batch(tasks, 1, "Sequential")

        # Run in parallel mode
        par_results, par_warnings, _ = run_parallel_batch(tasks, 2, "Parallel")

        # Same number of results
        assert len(seq_results) == len(par_results)
        assert len(seq_warnings) == len(par_warnings)

        # Same gene IDs processed
        seq_ids = {r.gene_id for r in seq_results}
        par_ids = {r.gene_id for r in par_results}
        assert seq_ids == par_ids

        # Same MK test values for each gene
        seq_by_id = {r.gene_id: r.result for r in seq_results}
        par_by_id = {r.gene_id: r.result for r in par_results}

        for gene_id in seq_ids:
            seq_mk = seq_by_id[gene_id]
            par_mk = par_by_id[gene_id]
            assert seq_mk.dn == par_mk.dn
            assert seq_mk.ds == par_mk.ds
            assert seq_mk.pn == par_mk.pn
            assert seq_mk.ps == par_mk.ps


class TestBatchTaskDataclass:
    """Tests for the BatchTask dataclass."""

    def test_default_values(self) -> None:
        """Test BatchTask has correct defaults."""
        task = BatchTask(file_path=Path("/test/file.fa"))

        assert task.ingroup_match is None
        assert task.outgroup_match is None
        assert task.outgroup_file is None
        assert task.reading_frame == 1
        assert task.use_asymptotic is False
        assert task.bins == 10
        assert task.bootstrap == 100
        assert task.pool_polymorphisms is False
        assert task.min_freq == 0.0
        assert task.extract_only is False

    def test_combined_mode_task(self) -> None:
        """Test creating a combined mode task."""
        task = BatchTask(
            file_path=Path("/test/alignment.fa"),
            ingroup_match="gamb",
            outgroup_match="afun",
            reading_frame=1,
        )

        assert task.ingroup_match == "gamb"
        assert task.outgroup_match == "afun"
        assert task.outgroup_file is None

    def test_separate_files_task(self) -> None:
        """Test creating a separate files mode task."""
        task = BatchTask(
            file_path=Path("/test/ingroup.fa"),
            outgroup_file=Path("/test/outgroup.fa"),
            reading_frame=1,
        )

        assert task.ingroup_match is None
        assert task.outgroup_file == Path("/test/outgroup.fa")


class TestWorkerResultDataclass:
    """Tests for the WorkerResult dataclass."""

    def test_success_result(self) -> None:
        """Test creating a successful result."""
        result = WorkerResult(gene_id="gene1", result={"test": "data"})

        assert result.gene_id == "gene1"
        assert result.result == {"test": "data"}
        assert result.error is None
        assert result.warning is None

    def test_error_result(self) -> None:
        """Test creating an error result."""
        result = WorkerResult(gene_id="gene1", error="Something went wrong")

        assert result.gene_id == "gene1"
        assert result.result is None
        assert result.error == "Something went wrong"

    def test_warning_result(self) -> None:
        """Test creating a warning result."""
        result = WorkerResult(gene_id="gene1", warning="Missing data")

        assert result.gene_id == "gene1"
        assert result.result is None
        assert result.warning == "Missing data"


class TestRunParallelBatchErrorFlag:
    """The empty-result exit code depends on whether the emptiness came from errors."""

    def test_reports_no_error_for_a_clean_run(self, tmp_path: Path) -> None:
        f = tmp_path / "gene1.fa"
        f.write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
        task = BatchTask(file_path=f, ingroup_match="speciesA", outgroup_match="speciesB")
        results, warnings, had_error = run_parallel_batch([task], 1, "Testing")
        assert len(results) == 1
        assert had_error is False

    def test_reports_no_error_when_input_is_only_unanalysable(self, tmp_path: Path) -> None:
        f = tmp_path / "gene1.fa"
        f.write_text(">speciesA_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
        task = BatchTask(file_path=f, ingroup_match="nomatch", outgroup_match="alsonomatch")
        results, warnings, had_error = run_parallel_batch([task], 1, "Testing")
        assert results == []
        assert warnings
        assert had_error is False

    def test_reports_an_error_when_a_gene_raises(self, tmp_path: Path) -> None:
        task = BatchTask(
            file_path=tmp_path / "missing.fa", ingroup_match="speciesA", outgroup_match="speciesB"
        )
        results, warnings, had_error = run_parallel_batch([task], 1, "Testing")
        assert results == []
        assert had_error is True
