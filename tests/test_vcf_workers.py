"""Tests for the VCF batch workers."""

from __future__ import annotations

import pickle
from pathlib import Path

import pytest

from mkado.analysis.asymptotic import AsymptoticMKResult, PolymorphismData
from mkado.analysis.imputed import ImputedMKResult
from mkado.analysis.mk_test import MKResult
from mkado.io.vcf import GeneStats
from mkado.vcf_workers import (
    VcfBatchChunk,
    VcfBatchTask,
    _process_single_gene,
    process_vcf_chunk,
    process_vcf_gene,
)

CHR1 = ["g_plus", "g_minus", "g_split", "g_phase1", "g_ncodon"]


def make_task(genome, spec_id: str, **overrides) -> VcfBatchTask:
    """Build a task for one synthetic gene, with any field overridden."""
    cds = genome.cds(spec_id)
    fields = dict(
        gene_id=cds.gene_id,
        transcript_id=cds.transcript_id,
        chrom=cds.chrom,
        exons=cds.exons,
        strand=cds.strand,
        phase=cds.phase,
        vcf_path=genome.ingroup_vcf,
        outgroup_vcf_path=genome.outgroup_vcf,
        ref_fasta_path=genome.ref_fasta,
    )
    return VcfBatchTask(**{**fields, **overrides})


def _counts(result: MKResult) -> tuple[int, int, int, int]:
    return (result.dn, result.ds, result.pn, result.ps)


def _key(wr):
    """The fields that must agree between the chunked and the per-gene worker."""
    counts = _counts(wr.result) if isinstance(wr.result, MKResult) else None
    p_value = wr.result.p_value if isinstance(wr.result, MKResult) else None
    return (wr.gene_id, wr.warning, wr.error, counts, p_value)


class TestVcfBatchTask:
    def test_defaults(self):
        task = VcfBatchTask(
            "g", "tx", "chr1", [(0, 9)], "+", 0, Path("in.vcf.gz"), None, Path("ref.fa")
        )
        assert task.code_table == 1
        assert task.min_freq == 0.0
        assert task.no_singletons is False
        assert task.use_asymptotic is False
        assert task.bins == 10
        assert task.bootstrap == 100
        assert task.use_imputed is False
        assert task.imputed_cutoff == 0.15
        assert task.extract_only is False
        assert task.ci_method == "monte-carlo"
        assert task.sfs_mode == "at"

    def test_picklable(self, genome):
        task = make_task(genome, "g_plus", use_asymptotic=True, bins=3)
        assert pickle.loads(pickle.dumps(task)) == task
        chunk = VcfBatchChunk(tasks=[task])
        assert pickle.loads(pickle.dumps(chunk)) == chunk


class TestProcessVcfGene:
    def test_standard_mk(self, genome):
        expected = genome.expected["g_plus"]
        wr = process_vcf_gene(make_task(genome, "g_plus"))
        assert wr.gene_id == "g_plus"
        assert wr.error is None
        assert isinstance(wr.result, MKResult)
        assert _counts(wr.result) == expected.counts
        assert wr.warning == expected.warning

    @pytest.mark.parametrize("gene_id", CHR1[1:])
    def test_standard_mk_without_skipped_sites(self, genome, gene_id):
        wr = process_vcf_gene(make_task(genome, gene_id))
        assert wr.warning is None
        assert _counts(wr.result) == genome.expected[gene_id].counts

    def test_extract_only_returns_polymorphism_data(self, genome):
        expected = genome.expected["g_plus"]
        wr = process_vcf_gene(make_task(genome, "g_plus", extract_only=True))
        assert isinstance(wr.result, PolymorphismData)
        assert wr.result.gene_id == "g_plus"
        assert (wr.result.dn, wr.result.ds) == (expected.dn, expected.ds)
        assert [kind for _, kind in wr.result.polymorphisms] == ["S", "N", "N"]
        assert wr.warning == expected.warning

    def test_extract_only_wins_over_other_modes(self, genome):
        task = make_task(genome, "g_plus", extract_only=True, use_asymptotic=True, use_imputed=True)
        assert isinstance(process_vcf_gene(task).result, PolymorphismData)

    @pytest.mark.parametrize("ci_method", ["monte-carlo", "bootstrap"])
    def test_asymptotic_per_gene(self, genome, ci_method):
        task = make_task(genome, "g_plus", use_asymptotic=True, bootstrap=2, ci_method=ci_method)
        wr = process_vcf_gene(task)
        assert wr.error is None
        result = wr.result
        assert isinstance(result, AsymptoticMKResult)
        assert (result.dn, result.ds, result.pn_total, result.ps_total) == (1, 2, 2, 1)
        assert result.num_genes == 1
        assert result.ci_method == ci_method
        assert result.alpha_asymptotic == pytest.approx(-3.0)

    def test_asymptotic_options_forwarded(self, genome):
        task = make_task(
            genome, "g_plus", use_asymptotic=True, bootstrap=1, bins=5, sfs_mode="above"
        )
        result = process_vcf_gene(task).result
        assert isinstance(result, AsymptoticMKResult)
        assert result.sfs_mode == "above"
        assert len(result.frequency_bins) <= 5

    def test_imputed_per_gene(self, genome):
        task = make_task(genome, "g_plus", use_imputed=True, imputed_cutoff=0.3, bootstrap=5)
        wr = process_vcf_gene(task)
        assert wr.error is None
        result = wr.result
        assert isinstance(result, ImputedMKResult)
        assert result.cutoff == 0.3
        assert (result.dn, result.ds, result.pn_total, result.ps_total) == (1, 2, 2, 1)

    def test_vertebrate_mito_code(self, genome):
        wr = process_vcf_gene(make_task(genome, "g_minus", code_table=2))
        assert _counts(wr.result) == (0, 2, 1, 1)

    def test_min_freq(self, genome):
        wr = process_vcf_gene(make_task(genome, "g_plus", min_freq=0.2))
        assert _counts(wr.result) == (1, 2, 2, 0)

    def test_no_singletons(self, genome):
        wr = process_vcf_gene(make_task(genome, "g_plus", no_singletons=True))
        assert _counts(wr.result) == (1, 2, 2, 0)

    def test_no_outgroup(self, genome):
        wr = process_vcf_gene(make_task(genome, "g_plus", outgroup_vcf_path=None))
        assert _counts(wr.result) == (0, 0, 2, 1)

    def test_bgzipped_reference(self, genome, bgzipped_ref):
        wr = process_vcf_gene(make_task(genome, "g_plus", ref_fasta_path=bgzipped_ref))
        assert wr.error is None
        assert _counts(wr.result) == genome.expected["g_plus"].counts

    def test_error_missing_reference(self, genome, tmp_path):
        wr = process_vcf_gene(make_task(genome, "g_plus", ref_fasta_path=tmp_path / "nope.fa"))
        assert wr.result is None
        assert wr.warning is None
        assert wr.error is not None
        assert wr.error.startswith("Error processing g_plus:")

    def test_error_chrom_absent_from_fasta(self, genome):
        wr = process_vcf_gene(make_task(genome, "g_nochrom"))
        assert wr.result is None
        assert wr.error is not None
        assert "Error processing g_nochrom:" in wr.error
        assert "chrZ" in wr.error


class TestWarningString:
    @pytest.mark.parametrize(
        ("indels", "multi", "expected"),
        [
            (3, 0, "g: skipped 3 indels"),
            (0, 2, "g: skipped 2 multi-allelic"),
            (3, 2, "g: skipped 3 indels, 2 multi-allelic"),
            (0, 0, None),
        ],
    )
    def test_warning_text(self, genome, monkeypatch, indels, multi, expected):
        def fake_extract(*args, **kwargs):
            data = PolymorphismData(polymorphisms=[], dn=0, ds=0, gene_id="g")
            return data, GeneStats(indels, multi, 0)

        monkeypatch.setattr("mkado.io.vcf.extract_gene_data", fake_extract)
        task = make_task(genome, "g_plus", gene_id="g", extract_only=True)
        wr = _process_single_gene(task, None, None, None)
        assert wr.error is None
        assert wr.warning == expected


class TestProcessVcfChunk:
    def test_chunk_equals_per_gene(self, genome):
        tasks = [make_task(genome, gene_id) for gene_id in CHR1]
        chunked = process_vcf_chunk(VcfBatchChunk(tasks=tasks))
        single = [process_vcf_gene(task) for task in tasks]
        assert [_key(wr) for wr in chunked] == [_key(wr) for wr in single]
        assert [wr.gene_id for wr in chunked] == CHR1

    def test_chunk_without_outgroup(self, genome):
        tasks = [make_task(genome, g, outgroup_vcf_path=None) for g in ("g_plus", "g_minus")]
        results = process_vcf_chunk(VcfBatchChunk(tasks=tasks))
        assert [_counts(wr.result) for wr in results] == [(0, 0, 2, 1), (0, 0, 1, 1)]

    def test_chunk_mixed_modes(self, genome):
        tasks = [
            make_task(genome, "g_split", extract_only=True),
            make_task(genome, "g_plus", use_asymptotic=True, bootstrap=1),
            make_task(genome, "g_minus", use_imputed=True, bootstrap=2),
            make_task(genome, "g_phase1"),
        ]
        results = process_vcf_chunk(VcfBatchChunk(tasks=tasks))
        assert [wr.error for wr in results] == [None] * 4
        assert [type(wr.result) for wr in results] == [
            PolymorphismData,
            AsymptoticMKResult,
            ImputedMKResult,
            MKResult,
        ]

    def test_chunk_error_isolated_to_gene(self, genome):
        tasks = [
            make_task(genome, "g_plus"),
            make_task(genome, "g_nochrom"),
            make_task(genome, "g_minus"),
        ]
        results = process_vcf_chunk(VcfBatchChunk(tasks=tasks))
        assert results[1].result is None
        assert "chrZ" in results[1].error
        assert _counts(results[0].result) == genome.expected["g_plus"].counts
        assert _counts(results[2].result) == genome.expected["g_minus"].counts

    def test_chunk_open_failure_propagates(self, genome, not_a_vcf):
        chunk = VcfBatchChunk(tasks=[make_task(genome, "g_plus", vcf_path=not_a_vcf)])
        with pytest.raises(OSError, match="not valid bcf or vcf"):
            process_vcf_chunk(chunk)
