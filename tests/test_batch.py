"""Tests for parallel batch processing."""

from pathlib import Path
from unittest.mock import patch

import pytest
from scipy.stats import false_discovery_control

from mkado.analysis.alpha_tg import AlphaTGResult
from mkado.analysis.asymptotic import AsymptoticMKResult, PolymorphismData
from mkado.analysis.imputed import ImputedMKResult
from mkado.analysis.mk_test import mk_test_from_counts
from mkado.analysis.polarized import PolarizedMKResult
from mkado.batch_workers import BatchTask, WorkerResult, process_gene
from mkado.cli import compute_adjusted_pvalues, get_worker_count, run_parallel_batch
from tests.builders import TwoGeneDirs


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

    def test_stop_blocked_pair_warning(self, tmp_path: Path) -> None:
        """Issue #90: a stop-blocked codon pair surfaces as a WorkerResult warning."""
        alignment = tmp_path / "gene1.fa"
        alignment.write_text(">gene1_speciesA_1\nAAA\n>gene1_speciesB_1\nTGG\n")

        task = BatchTask(
            file_path=alignment,
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            code_table=2,  # vertebrate mitochondrial
        )

        result = process_gene(task)

        assert result.error is None
        assert result.result is not None
        assert result.warning is not None
        assert "stop-free path" in result.warning

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

    def test_combined_mode_imputed(self, two_genes: TwoGeneDirs) -> None:
        """Combined-file imputed mode returns an ImputedMKResult."""
        task = BatchTask(
            file_path=two_genes.combined / "g1.fa",
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            use_imputed=True,
            bootstrap=10,
        )
        result = process_gene(task)
        assert result.error is None
        assert isinstance(result.result, ImputedMKResult)

    def test_combined_mode_polarized(self, two_genes: TwoGeneDirs) -> None:
        """Combined-file polarized mode reads the second outgroup by name pattern."""
        task = BatchTask(
            file_path=two_genes.combined / "g1.fa",
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            polarize_match="speciesC",
        )
        result = process_gene(task)
        assert result.error is None
        assert isinstance(result.result, PolarizedMKResult)
        assert (result.result.dn_outgroup, result.result.ds_outgroup) == (1, 1)

    def test_combined_mode_no_outgroup2(self, two_genes: TwoGeneDirs) -> None:
        """A polarize pattern that matches nothing is a warning, not an error."""
        task = BatchTask(
            file_path=two_genes.combined / "g1.fa",
            ingroup_match="speciesA",
            outgroup_match="speciesB",
            polarize_match="speciesZ",
        )
        result = process_gene(task)
        assert result.result is None
        assert result.warning == "No outgroup2 sequences in g1.fa"

    def test_separate_files_extract_only(self, two_genes: TwoGeneDirs) -> None:
        """Separate-files extract-only mode returns the gene's PolymorphismData."""
        task = BatchTask(
            file_path=two_genes.separate / "g1_ingroup.fa",
            outgroup_file=two_genes.separate / "g1_outgroup.fa",
            extract_only=True,
        )
        result = process_gene(task)
        assert result.error is None
        assert isinstance(result.result, PolymorphismData)
        assert (result.result.dn, result.result.ds) == (1, 1)

    def test_separate_files_asymptotic(self, two_genes: TwoGeneDirs) -> None:
        """Separate-files asymptotic mode returns an AsymptoticMKResult."""
        task = BatchTask(
            file_path=two_genes.separate / "g1_ingroup.fa",
            outgroup_file=two_genes.separate / "g1_outgroup.fa",
            use_asymptotic=True,
            bins=3,
            bootstrap=5,
        )
        result = process_gene(task)
        assert result.error is None
        assert isinstance(result.result, AsymptoticMKResult)

    def test_separate_files_imputed(self, two_genes: TwoGeneDirs) -> None:
        """Separate-files imputed mode returns an ImputedMKResult."""
        task = BatchTask(
            file_path=two_genes.separate / "g1_ingroup.fa",
            outgroup_file=two_genes.separate / "g1_outgroup.fa",
            use_imputed=True,
            bootstrap=10,
        )
        result = process_gene(task)
        assert result.error is None
        assert isinstance(result.result, ImputedMKResult)

    def test_separate_files_polarized(self, two_genes: TwoGeneDirs) -> None:
        """A second outgroup file switches separate-files mode to the polarized test."""
        task = BatchTask(
            file_path=two_genes.separate / "g1_ingroup.fa",
            outgroup_file=two_genes.separate / "g1_outgroup.fa",
            outgroup2_file=two_genes.separate / "g1_outgroup2.fa",
        )
        result = process_gene(task)
        assert result.error is None
        assert isinstance(result.result, PolarizedMKResult)
        assert (result.result.dn_outgroup, result.result.ds_outgroup) == (1, 1)


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

    def test_parallel_mode_collects_warnings_and_errors(self, tmp_path: Path) -> None:
        """A worker's warning and error reach the caller from the parallel path too."""
        for i in range(10):
            (tmp_path / f"gene{i}.fa").write_text(
                f">gene{i}_speciesA_1\nATGATGATG\n>gene{i}_speciesA_2\nATGCTGATG\n"
                f">gene{i}_speciesB_1\nATGGTGATG\n"
            )
        (tmp_path / "unmatched.fa").write_text(">other_1\nATGATGATG\n>speciesB_1\nATGGTGATG\n")
        files = [tmp_path / f"gene{i}.fa" for i in range(10)]
        files += [tmp_path / "unmatched.fa", tmp_path / "missing.fa"]
        tasks = [
            BatchTask(file_path=f, ingroup_match="speciesA", outgroup_match="speciesB")
            for f in files
        ]

        results, warnings, had_error = run_parallel_batch(tasks, 2, "Testing parallel")

        assert len(results) == 10
        assert had_error is True
        assert "Warning: No ingroup sequences in unmatched.fa" in warnings
        assert any(w.startswith("Error processing missing.fa") for w in warnings)

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


def _imputed_result(p_value: float) -> ImputedMKResult:
    return ImputedMKResult(
        alpha=0.5,
        p_value=p_value,
        pn_neutral=3.0,
        pwd=1.0,
        dn=10,
        ds=5,
        pn_total=4,
        ps_total=8,
        cutoff=0.15,
    )


def _asymptotic_result() -> AsymptoticMKResult:
    return AsymptoticMKResult(
        alpha_asymptotic=0.0,
        ci_low=0.0,
        ci_high=0.0,
        dn=1,
        ds=1,
    )


def _polarized_result(p_value_ingroup: float) -> PolarizedMKResult:
    return PolarizedMKResult(
        dn_ingroup=1,
        ds_ingroup=1,
        pn_ingroup=1,
        ps_ingroup=1,
        dn_outgroup=1,
        ds_outgroup=1,
        dn_unpolarized=1,
        ds_unpolarized=1,
        pn_unpolarized=1,
        ps_unpolarized=1,
        p_value_ingroup=p_value_ingroup,
        ni_ingroup=None,
        alpha_ingroup=None,
        dos_ingroup=None,
    )


class TestComputeAdjustedPvalues:
    """compute_adjusted_pvalues must not fake a p-value for types that lack one."""

    def test_imputed_results_get_real_bh_correction(self) -> None:
        results = [("geneA", _imputed_result(0.01)), ("geneB", _imputed_result(0.5))]
        adjusted = compute_adjusted_pvalues(results)
        assert adjusted != [1.0, 1.0]
        expected = list(false_discovery_control([0.01, 0.5], method="bh"))
        assert adjusted == expected

    def test_asymptotic_results_return_none(self) -> None:
        results = [("geneA", _asymptotic_result()), ("geneB", _asymptotic_result())]
        assert compute_adjusted_pvalues(results) is None

    def test_mk_results_still_adjusted(self) -> None:
        results = [
            ("geneA", mk_test_from_counts(dn=10, ds=5, pn=4, ps=8)),
            ("geneB", mk_test_from_counts(dn=1, ds=1, pn=1, ps=1)),
        ]
        adjusted = compute_adjusted_pvalues(results)
        expected_pvalues = [r.p_value for _, r in results]
        assert adjusted == list(false_discovery_control(expected_pvalues, method="bh"))

    def test_polarized_results_still_adjusted(self) -> None:
        results = [("geneA", _polarized_result(0.02)), ("geneB", _polarized_result(0.5))]
        adjusted = compute_adjusted_pvalues(results)
        expected = list(false_discovery_control([0.02, 0.5], method="bh"))
        assert adjusted == expected

    def test_unknown_result_type_raises(self) -> None:
        alpha_tg = AlphaTGResult(
            alpha_tg=0.0,
            ni_tg=0.0,
            ci_low=0.0,
            ci_high=0.0,
            num_genes=0,
            dn_total=0,
            ds_total=0,
            pn_total=0,
            ps_total=0,
        )
        with pytest.raises(TypeError, match="Unknown result type"):
            compute_adjusted_pvalues([("geneA", alpha_tg)])
