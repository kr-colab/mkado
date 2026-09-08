"""Tests for VCF input functionality."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pysam
import pytest

from mkado.core.cds import CdsRegion
from mkado.core.codons import GeneticCode
from mkado.io.vcf import (
    _complement_base,
    _open_vcf,
    _reconstruct_codon_with_sub,
    _ref_base_fetcher,
    extract_gene_data,
)

OUTGROUP = ["outgroup"]

needs_procfs = pytest.mark.skipif(not os.path.isdir("/proc/self/fd"), reason="needs procfs")


def _open_fd_count() -> int:
    """Descriptors open in this process. Raw pipe ends never raise ResourceWarning."""
    return len(os.listdir("/proc/self/fd"))


def _extract(genome, gene_id, **kwargs):
    """Run extract_gene_data on one gene of a dataset; keyword arguments override the call."""
    call = dict(
        vcf_path=genome.ingroup_vcf,
        outgroup_vcf_path=genome.outgroup_vcf,
        cds=genome.cds(gene_id),
        ref_fasta_path=genome.ref_fasta,
    )
    return extract_gene_data(**{**call, **kwargs})


def _assert_polys(actual, expected):
    assert [kind for _, kind in actual] == [kind for _, kind in expected]
    assert [freq for freq, _ in actual] == pytest.approx([freq for freq, _ in expected])


# ---- Fixtures for synthetic test data ----


@pytest.fixture
def synthetic_ref(tmp_path: Path) -> Path:
    """Create a tiny reference genome with one 'gene' (9bp = 3 codons).

    Sequence: ATG GCC AAA (Met-Ala-Lys)
    Positions: 0123456789 (0-based)
    """
    ref_path = tmp_path / "ref.fa"
    # chr1 is 30 bases to give room
    seq = "ATGGCCAAATTTTTTTTTTTTTTTTTTTTTT"
    ref_path.write_text(f">chr1\n{seq}\n")

    # Index it
    pysam.faidx(str(ref_path))

    return ref_path


@pytest.fixture
def simple_cds() -> CdsRegion:
    """CDS for a simple 3-codon gene on plus strand."""
    return CdsRegion(
        gene_id="test_gene",
        transcript_id="tx1",
        chrom="chr1",
        exons=[(0, 9)],  # ATG GCC AAA
        strand="+",
    )


@pytest.fixture
def ingroup_vcf_synonymous(tmp_path: Path, write_vcf) -> Path:
    """Ingroup VCF with one synonymous polymorphism.

    Site at position 5 (0-based): C->T in codon GCC -> GCT (both Ala).
    4 diploid samples: 3 hom-ref, 1 het => alt_freq = 1/8 = 0.125
    VCF is 1-based, so pos=6.
    """
    records = ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1"]
    return write_vcf(tmp_path / "ingroup_syn", records)


@pytest.fixture
def ingroup_vcf_nonsyn(tmp_path: Path, write_vcf) -> Path:
    """Ingroup VCF with one nonsynonymous polymorphism.

    Site at position 3 (0-based): G->A in codon GCC -> ACC (Ala -> Thr).
    4 diploid samples: 2 hom-ref, 2 het => alt_freq = 2/8 = 0.25
    VCF is 1-based, so pos=4.
    """
    records = ["chr1\t4\t.\tG\tA\t30\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1"]
    return write_vcf(tmp_path / "ingroup_nonsyn", records)


@pytest.fixture
def ingroup_vcf_codon2_het(tmp_path: Path, write_vcf) -> Path:
    """Ingroup VCF with one nonsynonymous polymorphism at codon 2, base 2.

    Site at position 4 (0-based): C->T in codon GCC -> GTC (Ala -> Val).
    4 diploid samples: 1 het, 3 hom-ref => alt_freq = 1/8 = 0.125
    VCF is 1-based, so pos=5.
    """
    records = ["chr1\t5\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0"]
    return write_vcf(tmp_path / "ingroup_codon2_het", records)


@pytest.fixture
def outgroup_vcf_divergent(tmp_path: Path, write_vcf) -> Path:
    """Outgroup VCF with a fixed difference at codon 3 (AAA -> AGA = Lys -> Arg).

    Site at position 7 (0-based): A->G. VCF pos=8.
    Single sample, hom-alt.
    """
    records = ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"]
    return write_vcf(tmp_path / "outgroup", records, samples=OUTGROUP)


@pytest.fixture
def outgroup_vcf_empty(tmp_path: Path, write_vcf) -> Path:
    """Outgroup VCF with no variants (all same as reference)."""
    return write_vcf(tmp_path / "outgroup_empty", [], samples=OUTGROUP)


@pytest.fixture
def ingroup_vcf_fixed_alt(tmp_path: Path, write_vcf) -> Path:
    """Ingroup fixed for G at position 7 (0-based): codon AAA -> AGA in every sample."""
    records = ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1"]
    return write_vcf(tmp_path / "ingroup_fixed_alt", records)


# htslib checks PL at header parse time and warns that it should be Number=G.
# The warning goes to file descriptor 2 during open, which is what _open_vcf
# captures. htslib prints it once per process, so this fixture must stay the
# only one in the suite whose header carries it.
PL_BAD_HEADER = '##FORMAT=<ID=PL,Number=1,Type=Integer,Description="Phred-scaled likelihoods">'


@pytest.fixture
def ingroup_vcf_pl_header(tmp_path: Path, write_vcf) -> Path:
    return write_vcf(tmp_path / "ingroup_pl", [], extra_header=[PL_BAD_HEADER])


# ---- Tests ----


class TestCodonReconstruction:
    def test_plus_strand_substitution(self, synthetic_ref, simple_cds):
        fasta = pysam.FastaFile(str(synthetic_ref))

        def ref_fetch(chrom, pos):
            return fasta.fetch(chrom, pos, pos + 1).upper()

        # Substitute position 5 (C->T) in codon GCC (index 1)
        ref_codon, alt_codon = _reconstruct_codon_with_sub(simple_cds, 1, ref_fetch, 5, "T")
        assert ref_codon == "GCC"
        assert alt_codon == "GCT"
        fasta.close()

    def test_nonsynonymous_substitution(self, synthetic_ref, simple_cds):
        fasta = pysam.FastaFile(str(synthetic_ref))

        def ref_fetch(chrom, pos):
            return fasta.fetch(chrom, pos, pos + 1).upper()

        # Substitute position 3 (G->A) in codon GCC (index 1)
        ref_codon, alt_codon = _reconstruct_codon_with_sub(simple_cds, 1, ref_fetch, 3, "A")
        assert ref_codon == "GCC"
        assert alt_codon == "ACC"
        fasta.close()


class TestExtractGeneData:
    def test_synonymous_polymorphism(
        self, synthetic_ref, simple_cds, ingroup_vcf_synonymous, outgroup_vcf_empty
    ):
        """One synonymous SNP should yield Ps=1, Pn=0."""
        poly_data, stats = extract_gene_data(
            vcf_path=ingroup_vcf_synonymous,
            outgroup_vcf_path=outgroup_vcf_empty,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )

        syn_count = sum(1 for _, t in poly_data.polymorphisms if t == "S")
        nonsyn_count = sum(1 for _, t in poly_data.polymorphisms if t == "N")
        assert syn_count == 1
        assert nonsyn_count == 0
        assert poly_data.dn == 0
        assert poly_data.ds == 0

    def test_nonsynonymous_polymorphism(
        self, synthetic_ref, simple_cds, ingroup_vcf_nonsyn, outgroup_vcf_empty
    ):
        """One nonsynonymous SNP should yield Pn=1, Ps=0."""
        poly_data, stats = extract_gene_data(
            vcf_path=ingroup_vcf_nonsyn,
            outgroup_vcf_path=outgroup_vcf_empty,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )

        syn_count = sum(1 for _, t in poly_data.polymorphisms if t == "S")
        nonsyn_count = sum(1 for _, t in poly_data.polymorphisms if t == "N")
        assert nonsyn_count == 1
        assert syn_count == 0

    def test_divergence_counting(
        self, synthetic_ref, simple_cds, ingroup_vcf_synonymous, outgroup_vcf_divergent
    ):
        """Outgroup difference at codon 3 (AAA->AGA, Lys->Arg) should give Dn=1."""
        poly_data, stats = extract_gene_data(
            vcf_path=ingroup_vcf_synonymous,
            outgroup_vcf_path=outgroup_vcf_divergent,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )

        # AAA -> AGA is nonsynonymous (Lys -> Arg)
        assert poly_data.dn == 1
        assert poly_data.ds == 0

    def test_no_outgroup(self, synthetic_ref, simple_cds, ingroup_vcf_synonymous):
        """Without outgroup, divergence should be 0."""
        poly_data, stats = extract_gene_data(
            vcf_path=ingroup_vcf_synonymous,
            outgroup_vcf_path=None,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )

        assert poly_data.dn == 0
        assert poly_data.ds == 0
        # Polymorphism should still be detected
        assert len(poly_data.polymorphisms) == 1

    def test_min_frequency_filter(
        self, synthetic_ref, simple_cds, ingroup_vcf_synonymous, outgroup_vcf_empty
    ):
        """With min_frequency > alt_freq, polymorphism should be filtered."""
        poly_data, stats = extract_gene_data(
            vcf_path=ingroup_vcf_synonymous,
            outgroup_vcf_path=outgroup_vcf_empty,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
            min_frequency=0.5,  # alt_freq is 0.125, should be filtered
        )
        assert len(poly_data.polymorphisms) == 0

    def test_gene_id_propagated(
        self, synthetic_ref, simple_cds, ingroup_vcf_synonymous, outgroup_vcf_empty
    ):
        poly_data, _ = extract_gene_data(
            vcf_path=ingroup_vcf_synonymous,
            outgroup_vcf_path=outgroup_vcf_empty,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )
        assert poly_data.gene_id == "test_gene"


class TestPolarization:
    def test_polarization_flips_frequency(self, synthetic_ref, simple_cds, tmp_path, write_vcf):
        """If outgroup carries ALT, derived freq should be 1 - alt_freq."""
        # Create ingroup with a SNP at pos 5 (C->T, synonymous)
        ingroup_records = ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1"]
        ingroup_vcf = write_vcf(tmp_path / "ig_polar", ingroup_records)

        # Outgroup also carries T at this position
        out_records = ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1"]
        outgroup_vcf = write_vcf(tmp_path / "og_polar", out_records, samples=OUTGROUP)

        poly_data, _ = extract_gene_data(
            vcf_path=ingroup_vcf,
            outgroup_vcf_path=outgroup_vcf,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )

        # alt_freq = 1/8 = 0.125, but outgroup has ALT so derived = 1 - 0.125 = 0.875
        assert len(poly_data.polymorphisms) == 1
        freq, _ = poly_data.polymorphisms[0]
        assert abs(freq - 0.875) < 0.01


# ---- Synthetic genome: parser branches ----


EXPECTED_CODONS = {
    "g_plus": "ATG GCC AAA TTC GGA CTG AGC TAC TAA",
    "g_minus": "ATG GAC AAG TTC",
    "g_split": "ATG AAA CCT GGG",
    "g_phase1": "ATG GTT CAC",
    "g_ncodon": "ATG CNT GGC",
    "tiny00": "ATG AAA TTT",
    "tiny11": "ATG AAA TTT",
}


class TestSyntheticGenome:
    def test_codons_match_design(self, genome):
        """Guard the reference literal: every gene must read as designed."""
        fasta = pysam.FastaFile(str(genome.ref_fasta))
        fetch = _ref_base_fetcher(fasta)
        for gene_id, expected in EXPECTED_CODONS.items():
            cds = genome.cds(gene_id)
            codons = [cds.extract_codon(i, fetch) for i in range(cds.num_codons())]
            assert " ".join(codons) == expected, gene_id
        fasta.close()

    @pytest.mark.parametrize(
        "gene_id",
        [
            pytest.param("g_plus", id="plus_strand_with_skipped_sites"),
            pytest.param("g_minus", id="minus_strand"),
            pytest.param("g_split", id="codon_spans_exon_junction"),
            pytest.param("g_phase1", id="phase_trimmed_base"),
            pytest.param("g_ncodon", id="ambiguous_reference_base"),
        ],
    )
    def test_counts_match_design(self, genome, gene_id):
        """Every record in the synthetic VCFs lands in the branch it was written for."""
        expected = genome.expected[gene_id]
        poly, stats = _extract(genome, gene_id)
        _assert_polys(poly.polymorphisms, expected.polymorphisms)
        assert (poly.dn, poly.ds) == (expected.dn, expected.ds)
        assert stats == expected.stats

    def test_g_minus_vertebrate_mito(self, genome):
        """Under table 2 ATA codes for Met, so the ATG to ATA difference becomes synonymous."""
        poly, _ = _extract(genome, "g_minus", genetic_code=GeneticCode(table_id=2))
        _assert_polys(poly.polymorphisms, genome.expected["g_minus"].polymorphisms)
        assert (poly.dn, poly.ds) == (0, 2)


class TestComplementBase:
    @pytest.mark.parametrize(
        ("base", "expected"), [("A", "T"), ("C", "G"), ("g", "c"), ("t", "a"), ("N", "N")]
    )
    def test_complement(self, base, expected):
        assert _complement_base(base) == expected


class TestHtslibWarningCapture:
    def test_bad_pl_header_logged(self, ingroup_vcf_pl_header, caplog):
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            vcf = _open_vcf(ingroup_vcf_pl_header)
        vcf.close()
        htslib = [r.message for r in caplog.records if r.message.startswith("htslib:")]
        assert htslib
        assert any("PL" in message for message in htslib)

    def test_clean_header_logs_nothing(self, genome, caplog):
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            vcf = _open_vcf(genome.ingroup_vcf)
        vcf.close()
        assert not [r for r in caplog.records if r.message.startswith("htslib:")]

    @needs_procfs
    def test_failed_open_closes_pipe(self, not_a_vcf):
        """The capture pipe must not leak a descriptor when htslib rejects the file."""
        before = _open_fd_count()
        with pytest.raises(OSError):
            _open_vcf(not_a_vcf)
        assert _open_fd_count() == before

    @needs_procfs
    def test_successful_open_closes_pipe(self, genome):
        before = _open_fd_count()
        _open_vcf(genome.ingroup_vcf).close()
        assert _open_fd_count() == before


class TestIngroupFiltering:
    def test_hom_alt_counts_two_alleles(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        records = ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1\t0/0\t0/0\t0/0"]
        ingroup = write_vcf(tmp_path / "hom_alt", records)
        poly, _ = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        _assert_polys(poly.polymorphisms, [(0.25, "S")])

    def test_fixed_alt_site_not_polymorphic(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        records = ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1"]
        ingroup = write_vcf(tmp_path / "fixed_alt", records)
        poly, _ = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []

    def test_shared_fixed_alt_is_not_divergence(
        self, synthetic_ref, simple_cds, ingroup_vcf_fixed_alt, outgroup_vcf_divergent
    ):
        """Ingroup and outgroup both carry G at codon 3; there is no difference between them."""
        poly, _ = extract_gene_data(
            ingroup_vcf_fixed_alt, outgroup_vcf_divergent, simple_cds, synthetic_ref
        )
        assert poly.polymorphisms == []
        assert (poly.dn, poly.ds) == (0, 0)

    def test_fixed_alt_against_reference_outgroup_is_divergence(
        self, synthetic_ref, simple_cds, ingroup_vcf_fixed_alt, outgroup_vcf_empty
    ):
        """The ingroup is fixed for G at codon 3 and the outgroup matches the reference A."""
        poly, _ = extract_gene_data(
            ingroup_vcf_fixed_alt, outgroup_vcf_empty, simple_cds, synthetic_ref
        )
        assert poly.polymorphisms == []
        assert (poly.dn, poly.ds) == (1, 0)


class TestOutgroupParsing:
    def test_sites_only_outgroup_ignored(self, genome, outgroup_vcf_sites_only):
        """With no outgroup alleles, only the ingroup's fixed ALT at codon 3 is a difference."""
        poly, _ = _extract(genome, "g_plus", outgroup_vcf_path=outgroup_vcf_sites_only)
        assert (poly.dn, poly.ds) == (1, 0)
        # Without an outgroup allele at codon 5 the polymorphism there is not polarized.
        _assert_polys(poly.polymorphisms, [(0.125, "S"), (0.25, "N"), (0.25, "N")])

    def test_het_outgroup_not_counted_as_alt(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        """A heterozygous outgroup call is unresolved, like a missing record, not ALT."""
        ingroup = write_vcf(tmp_path / "het_in", [])
        outgroup = write_vcf(
            tmp_path / "het_out", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"], samples=OUTGROUP
        )
        poly, _ = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)

    def test_het_outgroup_does_not_polarize(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, ingroup_vcf_codon2_het
    ):
        """A heterozygous outgroup call cannot resolve the ancestral allele, so the
        ingroup polymorphism keeps its raw ALT frequency instead of being flipped."""
        outgroup = write_vcf(
            tmp_path / "poly_out", ["chr1\t5\t.\tC\tT\t30\tPASS\t.\tGT\t0/1"], samples=OUTGROUP
        )
        poly, _ = extract_gene_data(ingroup_vcf_codon2_het, outgroup, simple_cds, synthetic_ref)
        _assert_polys(poly.polymorphisms, [(0.125, "N")])

    def test_indel_multiallelic_symbolic_ignored(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        ingroup = write_vcf(tmp_path / "skip_in", [])
        records = [
            "chr1\t2\t.\tTG\tT\t30\tPASS\t.\tGT\t1/1",
            "chr1\t4\t.\tG\tA,T\t30\tPASS\t.\tGT\t1/1",
            "chr1\t7\t.\tA\t<DEL>\t30\tPASS\tSVTYPE=DEL\tGT\t1/1",
        ]
        outgroup = write_vcf(tmp_path / "skip_out", records, samples=OUTGROUP)
        poly, _ = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)

    def test_ref_mismatch_identical_codon_skipped(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """The record's ALT equals the FASTA base, so the reconstructed codons are identical."""
        ingroup = write_vcf(tmp_path / "same_in", [])
        outgroup = write_vcf(
            tmp_path / "same_out", ["chr1\t6\t.\tT\tC\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, _ = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)

    def test_third_allele_at_polymorphic_site_is_ignored(
        self, synthetic_ref, simple_cds, ingroup_vcf_nonsyn, write_vcf, tmp_path
    ):
        """The ingroup is G/A at position 3 and the outgroup carries T: no fixed difference."""
        outgroup = write_vcf(
            tmp_path / "third_out", ["chr1\t4\t.\tG\tT\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, _ = extract_gene_data(ingroup_vcf_nonsyn, outgroup, simple_cds, synthetic_ref)
        _assert_polys(poly.polymorphisms, [(0.25, "N")])
        assert (poly.dn, poly.ds) == (0, 0)

    def test_created_stop_codon_skipped(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        """AAA to TAA would create a stop, so the difference is not counted."""
        ingroup = write_vcf(tmp_path / "stop_in", [])
        outgroup = write_vcf(
            tmp_path / "stop_out", ["chr1\t7\t.\tA\tT\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, _ = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)


class TestPartiallyPolymorphicCodon:
    def test_whole_codon_discarded_when_ingroup_polymorphic(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, ingroup_vcf_codon2_het
    ):
        """Codon GCC: the ingroup is polymorphic at base 2 and the outgroup differs at bases 2 and 3.

        A polymorphic base voids the whole codon as a fixed difference, matching the
        FASTA path's requirement of one clean codon per group, so base 3's clean
        divergence does not count either. The polymorphism at base 2 is still polarized,
        since polarization reads the outgroup's raw genotype, not the divergence codon.
        """
        outgroup = write_vcf(
            tmp_path / "part_out",
            [
                "chr1\t5\t.\tC\tT\t30\tPASS\t.\tGT\t1/1",
                "chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1",
            ],
            samples=OUTGROUP,
        )
        poly, _ = extract_gene_data(ingroup_vcf_codon2_het, outgroup, simple_cds, synthetic_ref)
        _assert_polys(poly.polymorphisms, [(0.875, "N")])
        assert (poly.dn, poly.ds) == (0, 0)


class TestHandles:
    def test_preopened_handles_match_paths(self, genome):
        ingroup = _open_vcf(genome.ingroup_vcf)
        outgroup = _open_vcf(genome.outgroup_vcf)
        fasta = pysam.FastaFile(str(genome.ref_fasta))
        try:
            via_handles, stats_handles = extract_gene_data(
                vcf_path=Path("unused.vcf.gz"),
                outgroup_vcf_path=Path("unused_out.vcf.gz"),
                cds=genome.cds("g_plus"),
                ref_fasta_path=Path("unused.fa"),
                ingroup_vcf=ingroup,
                outgroup_vcf=outgroup,
                ref_fasta=fasta,
            )
            # The caller owns the handles, so they stay open after the call.
            assert fasta.fetch("chr1", 0, 3) == "ATG"
        finally:
            ingroup.close()
            outgroup.close()
            fasta.close()

        via_paths, stats_paths = _extract(genome, "g_plus")
        assert via_handles == via_paths
        assert stats_handles == stats_paths


class TestQueryFailures:
    def test_unindexed_ingroup_logged_and_empty(self, genome, ingroup_vcf_plain, caplog):
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            poly, _ = _extract(genome, "g_plus", vcf_path=ingroup_vcf_plain)
        assert poly.polymorphisms == []
        assert any("ingroup VCF query failed for g_plus" in r.message for r in caplog.records)

    def test_unindexed_outgroup_logged_and_treated_as_reference(
        self, genome, outgroup_vcf_plain, caplog
    ):
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            poly, _ = _extract(genome, "g_plus", outgroup_vcf_path=outgroup_vcf_plain)
        assert (poly.dn, poly.ds) == (1, 0)
        assert any("outgroup VCF query failed for g_plus" in r.message for r in caplog.records)


class TestSingletons:
    def test_no_singletons_removes_singletons(self, genome):
        """Four diploid samples make 0.125 the singleton frequency."""
        poly, _ = _extract(genome, "g_plus", no_singletons=True)
        _assert_polys(poly.polymorphisms, [(0.25, "N"), (0.75, "N")])

    def test_no_singletons_keeps_higher_min_frequency(self, genome):
        poly, _ = _extract(genome, "g_plus", no_singletons=True, min_frequency=0.5)
        _assert_polys(poly.polymorphisms, [(0.75, "N")])

    def test_no_singletons_removes_polarized_singleton(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """Seven of eight ALT with the outgroup on ALT makes REF the derived singleton."""
        ingroup = write_vcf(
            tmp_path / "pol_in", ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1\t1/1\t1/1\t0/1"]
        )
        outgroup = write_vcf(
            tmp_path / "pol_out", ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        kept, _ = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        _assert_polys(kept.polymorphisms, [(0.125, "S")])
        dropped, _ = extract_gene_data(
            ingroup, outgroup, simple_cds, synthetic_ref, no_singletons=True
        )
        assert dropped.polymorphisms == []

    def test_no_singletons_with_missing_genotypes(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """A single ALT copy is a singleton whatever the number of called samples."""
        ingroup = write_vcf(
            tmp_path / "miss_in", ["chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t./.\t0/1\t0/0\t0/0"]
        )
        poly, _ = extract_gene_data(ingroup, None, simple_cds, synthetic_ref, no_singletons=True)
        assert poly.polymorphisms == []

    def test_no_singletons_without_snps(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        ingroup = write_vcf(tmp_path / "empty_in", [])
        poly, _ = extract_gene_data(ingroup, None, simple_cds, synthetic_ref, no_singletons=True)
        assert poly.polymorphisms == []


class TestFrequencyFilter:
    def test_min_frequency_on_g_plus(self, genome):
        poly, _ = _extract(genome, "g_plus", min_frequency=0.2)
        _assert_polys(poly.polymorphisms, [(0.25, "N"), (0.75, "N")])

    def test_site_at_min_frequency_is_kept(self, genome):
        """Only sites below the minimum are dropped, as the option's documentation says."""
        poly, _ = _extract(genome, "g_plus", min_frequency=0.125)
        _assert_polys(poly.polymorphisms, genome.expected["g_plus"].polymorphisms)
