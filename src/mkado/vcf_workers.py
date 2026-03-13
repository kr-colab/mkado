"""Worker functions for parallel VCF-based batch processing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mkado.batch_workers import WorkerResult


@dataclass
class VcfBatchTask:
    """Task specification for VCF-based gene processing.

    Designed to be picklable for multiprocessing.
    """

    gene_id: str
    transcript_id: str
    chrom: str
    exons: list[tuple[int, int]]
    strand: str
    phase: int

    vcf_path: Path
    outgroup_vcf_path: Path | None
    ref_fasta_path: Path

    code_table: int = 1
    min_freq: float = 0.0
    no_singletons: bool = False

    # Analysis flags
    use_asymptotic: bool = False
    bins: int = 10
    bootstrap: int = 100
    use_imputed: bool = False
    imputed_cutoff: float = 0.15
    extract_only: bool = False


@dataclass
class VcfBatchChunk:
    """A chunk of genes to process in a single worker invocation.

    All tasks in a chunk must share the same file paths so that
    handles can be opened once and reused.
    """

    tasks: list[VcfBatchTask]


def _process_single_gene(
    task: VcfBatchTask,
    ingroup_vcf: object,
    outgroup_vcf: object | None,
    ref_fasta: object,
) -> WorkerResult:
    """Process a single gene using pre-opened file handles."""
    from mkado.core.cds import CdsRegion
    from mkado.core.codons import GeneticCode
    from mkado.io.vcf import extract_gene_data

    try:
        cds = CdsRegion(
            gene_id=task.gene_id,
            transcript_id=task.transcript_id,
            chrom=task.chrom,
            exons=task.exons,
            strand=task.strand,
            phase=task.phase,
        )

        genetic_code = GeneticCode(table_id=task.code_table) if task.code_table != 1 else None

        poly_data, stats = extract_gene_data(
            vcf_path=task.vcf_path,
            outgroup_vcf_path=task.outgroup_vcf_path,
            cds=cds,
            ref_fasta_path=task.ref_fasta_path,
            genetic_code=genetic_code,
            min_frequency=task.min_freq,
            no_singletons=task.no_singletons,
            ingroup_vcf=ingroup_vcf,
            outgroup_vcf=outgroup_vcf,
            ref_fasta=ref_fasta,
        )

        # Build warning if many sites skipped
        warning = None
        skipped = stats.skipped_indels + stats.skipped_multiallelic
        if skipped > 0:
            parts = []
            if stats.skipped_indels > 0:
                parts.append(f"{stats.skipped_indels} indels")
            if stats.skipped_multiallelic > 0:
                parts.append(f"{stats.skipped_multiallelic} multi-allelic")
            warning = f"{task.gene_id}: skipped {', '.join(parts)}"

        # If extract_only, return PolymorphismData directly
        if task.extract_only:
            return WorkerResult(
                gene_id=task.gene_id,
                result=poly_data,
                warning=warning,
            )

        # Per-gene asymptotic
        if task.use_asymptotic:
            from mkado.analysis.asymptotic import asymptotic_mk_test_aggregated

            result = asymptotic_mk_test_aggregated(
                gene_data=[poly_data],
                num_bins=task.bins,
                ci_replicates=task.bootstrap * 100,
            )
            return WorkerResult(gene_id=task.gene_id, result=result, warning=warning)

        # Per-gene imputed
        if task.use_imputed:
            from mkado.analysis.imputed import imputed_mk_test

            result = imputed_mk_test(poly_data, cutoff=task.imputed_cutoff)
            return WorkerResult(gene_id=task.gene_id, result=result, warning=warning)

        # Standard MK from counts
        from mkado.analysis.mk_test import mk_test_from_counts

        pn = sum(1 for _, t in poly_data.polymorphisms if t == "N")
        ps = sum(1 for _, t in poly_data.polymorphisms if t == "S")

        result = mk_test_from_counts(
            dn=poly_data.dn,
            ds=poly_data.ds,
            pn=pn,
            ps=ps,
        )
        return WorkerResult(gene_id=task.gene_id, result=result, warning=warning)

    except Exception as e:
        return WorkerResult(
            gene_id=task.gene_id,
            error=f"Error processing {task.gene_id}: {e}",
        )


def process_vcf_chunk(chunk: VcfBatchChunk) -> list[WorkerResult]:
    """Process multiple genes with shared file handles.

    Opens VCF and FASTA handles once, then processes all genes in the chunk.
    This avoids the per-gene overhead of opening file handles.
    """
    import pysam

    from mkado.io.vcf import _open_vcf

    first = chunk.tasks[0]
    ingroup_vcf = _open_vcf(first.vcf_path)
    outgroup_vcf = _open_vcf(first.outgroup_vcf_path) if first.outgroup_vcf_path else None
    ref_fasta = pysam.FastaFile(str(first.ref_fasta_path))

    results = []
    for task in chunk.tasks:
        result = _process_single_gene(task, ingroup_vcf, outgroup_vcf, ref_fasta)
        results.append(result)

    ingroup_vcf.close()
    if outgroup_vcf:
        outgroup_vcf.close()
    ref_fasta.close()

    return results


def process_vcf_gene(task: VcfBatchTask) -> WorkerResult:
    """Process a single gene from VCF data.

    This function is designed to be called from a ProcessPoolExecutor.

    Args:
        task: VcfBatchTask with all parameters.

    Returns:
        WorkerResult with the analysis result or error message.
    """
    from mkado.core.cds import CdsRegion
    from mkado.core.codons import GeneticCode
    from mkado.io.vcf import extract_gene_data

    try:
        # Reconstruct CdsRegion from task fields
        cds = CdsRegion(
            gene_id=task.gene_id,
            transcript_id=task.transcript_id,
            chrom=task.chrom,
            exons=task.exons,
            strand=task.strand,
            phase=task.phase,
        )

        genetic_code = GeneticCode(table_id=task.code_table) if task.code_table != 1 else None

        poly_data, stats = extract_gene_data(
            vcf_path=task.vcf_path,
            outgroup_vcf_path=task.outgroup_vcf_path,
            cds=cds,
            ref_fasta_path=task.ref_fasta_path,
            genetic_code=genetic_code,
            min_frequency=task.min_freq,
            no_singletons=task.no_singletons,
        )

        # Build warning if many sites skipped
        warning = None
        skipped = stats.skipped_indels + stats.skipped_multiallelic
        if skipped > 0:
            parts = []
            if stats.skipped_indels > 0:
                parts.append(f"{stats.skipped_indels} indels")
            if stats.skipped_multiallelic > 0:
                parts.append(f"{stats.skipped_multiallelic} multi-allelic")
            warning = f"{task.gene_id}: skipped {', '.join(parts)}"

        # If extract_only, return PolymorphismData directly
        if task.extract_only:
            return WorkerResult(
                gene_id=task.gene_id,
                result=poly_data,
                warning=warning,
            )

        # Per-gene asymptotic
        if task.use_asymptotic:
            from mkado.analysis.asymptotic import asymptotic_mk_test_aggregated

            result = asymptotic_mk_test_aggregated(
                gene_data=[poly_data],
                num_bins=task.bins,
                ci_replicates=task.bootstrap * 100,
            )
            return WorkerResult(gene_id=task.gene_id, result=result, warning=warning)

        # Per-gene imputed
        if task.use_imputed:
            from mkado.analysis.imputed import imputed_mk_test

            result = imputed_mk_test(poly_data, cutoff=task.imputed_cutoff)
            return WorkerResult(gene_id=task.gene_id, result=result, warning=warning)

        # Standard MK from counts
        from mkado.analysis.mk_test import mk_test_from_counts

        pn = sum(1 for _, t in poly_data.polymorphisms if t == "N")
        ps = sum(1 for _, t in poly_data.polymorphisms if t == "S")

        result = mk_test_from_counts(
            dn=poly_data.dn,
            ds=poly_data.ds,
            pn=pn,
            ps=ps,
        )
        return WorkerResult(gene_id=task.gene_id, result=result, warning=warning)

    except Exception as e:
        return WorkerResult(
            gene_id=task.gene_id,
            error=f"Error processing {task.gene_id}: {e}",
        )
