"""Worker functions for parallel batch processing.

This module contains pure, picklable functions for use with ProcessPoolExecutor.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class BatchTask:
    """Task specification for batch processing.

    This dataclass holds all information needed to process a single gene,
    and is designed to be picklable for multiprocessing.
    """

    file_path: Path
    """Path to the alignment file (combined mode) or ingroup file (separate mode)."""

    # Combined file mode settings
    ingroup_match: str | None = None
    """Pattern to filter ingroup sequences (combined file mode)."""

    outgroup_match: str | None = None
    """Pattern to filter outgroup sequences (combined file mode)."""

    polarize_match: str | None = None
    """Pattern to filter second outgroup sequences (combined file polarized mode)."""

    # Separate files mode settings
    outgroup_file: Path | None = None
    """Path to outgroup file (separate files mode)."""

    outgroup2_file: Path | None = None
    """Path to second outgroup file (separate files polarized mode)."""

    # Analysis settings
    reading_frame: int = 1
    """Reading frame (1, 2, or 3)."""

    use_asymptotic: bool = False
    """Use asymptotic MK test instead of standard."""

    bins: int = 10
    """Number of frequency bins (asymptotic only)."""

    bootstrap: int = 100
    """Number of bootstrap replicates (asymptotic only)."""

    pool_polymorphisms: bool = False
    """Pool polymorphisms from both populations."""

    min_freq: float = 0.0
    """Minimum derived allele frequency for polymorphisms."""

    no_singletons: bool = False
    """Exclude singletons (automatically calculate 1/n threshold)."""

    use_imputed: bool = False
    """Use imputed MK test instead of standard."""

    imputed_cutoff: float = 0.15
    """DAF cutoff for imputed MK test."""

    extract_only: bool = False
    """Only extract polymorphism data (for aggregated asymptotic mode)."""

    code_table: int = 1
    """NCBI genetic code table ID."""


@dataclass
class WorkerResult:
    """Result from a worker function.

    Designed to be picklable for return from worker processes.
    """

    gene_id: str
    """Identifier for the gene (usually filename stem)."""

    result: Any | None = None
    """The analysis result (MKResult, PolarizedMKResult, AsymptoticMKResult, or PolymorphismData)."""

    error: str | None = None
    """Error message if processing failed."""

    warning: str | None = None
    """Warning message if there were issues."""


def process_gene(task: BatchTask) -> WorkerResult:
    """Process a single gene - pure function for multiprocessing.

    This function is designed to be called from a ProcessPoolExecutor.
    It contains no side effects and returns a picklable result.

    Args:
        task: BatchTask containing all parameters for processing

    Returns:
        WorkerResult with the analysis result or error message
    """
    # Import here to avoid pickling issues
    from mkado.analysis.asymptotic import (
        asymptotic_mk_test,
        extract_polymorphism_data,
    )
    from mkado.analysis.imputed import imputed_mk_test
    from mkado.analysis.mk_test import mk_test
    from mkado.analysis.polarized import polarized_mk_test
    from mkado.core.codons import GeneticCode
    from mkado.core.sequences import SequenceSet

    gene_id = task.file_path.stem
    genetic_code = GeneticCode(table_id=task.code_table) if task.code_table != 1 else None

    try:
        # Determine mode: combined file or separate files
        if task.ingroup_match is not None:
            # Combined file mode
            all_seqs = SequenceSet.from_fasta(
                task.file_path, reading_frame=task.reading_frame
            )
            ingroup_seqs = all_seqs.filter_by_name(task.ingroup_match)
            outgroup_seqs = all_seqs.filter_by_name(task.outgroup_match)

            if len(ingroup_seqs) == 0:
                return WorkerResult(
                    gene_id=gene_id,
                    warning=f"No ingroup sequences in {task.file_path.name}",
                )
            if len(outgroup_seqs) == 0:
                return WorkerResult(
                    gene_id=gene_id,
                    warning=f"No outgroup sequences in {task.file_path.name}",
                )

            # Calculate singleton threshold if --no-singletons
            min_freq = task.min_freq
            if task.no_singletons:
                n_samples = len(ingroup_seqs)
                if task.pool_polymorphisms:
                    n_samples += len(outgroup_seqs)
                min_freq = 1.0 / n_samples

            # Extract only mode (for aggregated asymptotic or alpha-tg)
            if task.extract_only:
                result = extract_polymorphism_data(
                    ingroup=ingroup_seqs,
                    outgroup=outgroup_seqs,
                    reading_frame=task.reading_frame,
                    pool_polymorphisms=task.pool_polymorphisms,
                    gene_id=gene_id,
                    min_frequency=min_freq,
                    genetic_code=genetic_code,
                )
                return WorkerResult(gene_id=gene_id, result=result)

            # Asymptotic mode
            if task.use_asymptotic:
                result = asymptotic_mk_test(
                    ingroup=ingroup_seqs,
                    outgroup=outgroup_seqs,
                    reading_frame=task.reading_frame,
                    num_bins=task.bins,
                    bootstrap_replicates=task.bootstrap,
                    pool_polymorphisms=task.pool_polymorphisms,
                    genetic_code=genetic_code,
                )
                return WorkerResult(gene_id=gene_id, result=result)

            # Imputed mode
            if task.use_imputed:
                poly_data = extract_polymorphism_data(
                    ingroup=ingroup_seqs,
                    outgroup=outgroup_seqs,
                    reading_frame=task.reading_frame,
                    pool_polymorphisms=task.pool_polymorphisms,
                    gene_id=gene_id,
                    genetic_code=genetic_code,
                )
                result = imputed_mk_test(poly_data, cutoff=task.imputed_cutoff)
                return WorkerResult(gene_id=gene_id, result=result)

            # Polarized mode
            if task.polarize_match:
                outgroup2_seqs = all_seqs.filter_by_name(task.polarize_match)
                if len(outgroup2_seqs) == 0:
                    return WorkerResult(
                        gene_id=gene_id,
                        warning=f"No outgroup2 sequences in {task.file_path.name}",
                    )
                result = polarized_mk_test(
                    ingroup=ingroup_seqs,
                    outgroup1=outgroup_seqs,
                    outgroup2=outgroup2_seqs,
                    reading_frame=task.reading_frame,
                    pool_polymorphisms=task.pool_polymorphisms,
                    min_frequency=min_freq,
                    genetic_code=genetic_code,
                )
                return WorkerResult(gene_id=gene_id, result=result)

            # Standard MK test
            result = mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=task.reading_frame,
                pool_polymorphisms=task.pool_polymorphisms,
                min_frequency=min_freq,
                genetic_code=genetic_code,
            )
            return WorkerResult(gene_id=gene_id, result=result)

        else:
            # Separate files mode
            if task.outgroup_file is None:
                return WorkerResult(
                    gene_id=gene_id,
                    warning=f"No outgroup found for {task.file_path.name}",
                )

            # Calculate singleton threshold if --no-singletons
            min_freq = task.min_freq
            if task.no_singletons:
                ingroup_seqs = SequenceSet.from_fasta(
                    task.file_path, reading_frame=task.reading_frame
                )
                n_samples = len(ingroup_seqs)
                if task.pool_polymorphisms:
                    outgroup_seqs = SequenceSet.from_fasta(
                        task.outgroup_file, reading_frame=task.reading_frame
                    )
                    n_samples += len(outgroup_seqs)
                min_freq = 1.0 / n_samples

            # Extract only mode (for aggregated asymptotic or alpha-tg)
            if task.extract_only:
                result = extract_polymorphism_data(
                    ingroup=task.file_path,
                    outgroup=task.outgroup_file,
                    reading_frame=task.reading_frame,
                    pool_polymorphisms=task.pool_polymorphisms,
                    gene_id=gene_id,
                    min_frequency=min_freq,
                    genetic_code=genetic_code,
                )
                return WorkerResult(gene_id=gene_id, result=result)

            # Asymptotic mode
            if task.use_asymptotic:
                result = asymptotic_mk_test(
                    ingroup=task.file_path,
                    outgroup=task.outgroup_file,
                    reading_frame=task.reading_frame,
                    num_bins=task.bins,
                    bootstrap_replicates=task.bootstrap,
                    pool_polymorphisms=task.pool_polymorphisms,
                    genetic_code=genetic_code,
                )
                return WorkerResult(gene_id=gene_id, result=result)

            # Imputed mode
            if task.use_imputed:
                poly_data = extract_polymorphism_data(
                    ingroup=task.file_path,
                    outgroup=task.outgroup_file,
                    reading_frame=task.reading_frame,
                    pool_polymorphisms=task.pool_polymorphisms,
                    gene_id=gene_id,
                    genetic_code=genetic_code,
                )
                result = imputed_mk_test(poly_data, cutoff=task.imputed_cutoff)
                return WorkerResult(gene_id=gene_id, result=result)

            # Polarized mode
            if task.outgroup2_file is not None:
                result = polarized_mk_test(
                    ingroup=task.file_path,
                    outgroup1=task.outgroup_file,
                    outgroup2=task.outgroup2_file,
                    reading_frame=task.reading_frame,
                    pool_polymorphisms=task.pool_polymorphisms,
                    min_frequency=min_freq,
                    genetic_code=genetic_code,
                )
                return WorkerResult(gene_id=gene_id, result=result)

            # Standard MK test
            result = mk_test(
                ingroup=task.file_path,
                outgroup=task.outgroup_file,
                reading_frame=task.reading_frame,
                pool_polymorphisms=task.pool_polymorphisms,
                min_frequency=min_freq,
                genetic_code=genetic_code,
            )
            return WorkerResult(gene_id=gene_id, result=result)

    except Exception as e:
        return WorkerResult(
            gene_id=gene_id,
            error=f"Error processing {task.file_path.name}: {e}",
        )
