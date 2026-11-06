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
    GeneStats,
    VcfQueryError,
    _complement_base,
    _htslib_stderr_to_log,
    _open_vcf,
    _reconstruct_codon_with_sub,
    _ref_base_fetcher,
    extract_gene_data,
    genotyped_samples,
    probe_vcf,
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
    def test_failed_open_closes_the_capture(self, not_a_vcf):
        """The capture must not leak a descriptor when htslib rejects the file."""
        before = _open_fd_count()
        with pytest.raises(OSError):
            _open_vcf(not_a_vcf)
        assert _open_fd_count() == before

    @needs_procfs
    def test_successful_open_closes_the_capture(self, genome):
        before = _open_fd_count()
        _open_vcf(genome.ingroup_vcf).close()
        assert _open_fd_count() == before

    def test_query_warning_arises_inside_the_query(
        self, synthetic_ref, simple_cds, ingroup_vcf_xq, caplog
    ):
        """The undeclared tag is reported when a record is parsed, inside the region
        query, not at open. That is why the workers hold the capture around the
        extraction rather than around the open alone."""
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            vcf = _open_vcf(ingroup_vcf_xq)
            assert not [r for r in caplog.records if "XQ" in r.message]
            with _htslib_stderr_to_log():
                extract_gene_data(ingroup_vcf_xq, None, simple_cds, synthetic_ref, ingroup_vcf=vcf)
            vcf.close()
        htslib = [r.message for r in caplog.records if r.message.startswith("htslib:")]
        assert any("XQ" in message for message in htslib)

    def test_capture_has_no_size_ceiling(self, caplog):
        """A hundred kilobytes of htslib lines written inside the capture must all reach the logger."""
        line = b"[W::test] " + b"x" * 89 + b"\n"
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            with _htslib_stderr_to_log():
                for _ in range(1000):
                    os.write(2, line)
        assert sum(1 for r in caplog.records if r.message.startswith("htslib:")) == 1000

    def test_other_output_in_the_window_is_passed_back_through(self, capfd, caplog):
        """A progress bar redraw that lands in the window goes back to the terminal, not the log."""
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            with _htslib_stderr_to_log():
                os.write(2, b"\x1b[2K progress 40%\n")
                os.write(2, b"[W::vcf_parse] an htslib line\n")
        assert capfd.readouterr().err == "\x1b[2K progress 40%\n"
        assert [r.message for r in caplog.records] == ["htslib: [W::vcf_parse] an htslib line"]

    def test_htslib_line_after_an_unterminated_redraw_is_still_logged(self, capfd, caplog):
        """A redraw with no newline does not swallow the htslib line written right after it."""
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            with _htslib_stderr_to_log():
                os.write(2, b"\x1b[2K progress 40%")
                os.write(2, b"[W::vcf_parse] an htslib line\n")
        assert capfd.readouterr().err == "\x1b[2K progress 40%"
        assert [r.message for r in caplog.records] == ["htslib: [W::vcf_parse] an htslib line"]

    def test_nested_capture_drains_once_at_the_outer_level(self, caplog):
        """An open inside a captured extraction must not start a second capture."""
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            with _htslib_stderr_to_log():
                with _htslib_stderr_to_log():
                    os.write(2, b"[W::inner] one\n")
                assert not caplog.records
                os.write(2, b"[W::outer] two\n")
        assert [r.message for r in caplog.records] == [
            "htslib: [W::inner] one",
            "htslib: [W::outer] two",
        ]


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

    def test_stop_blocked_pair_skipped_and_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """AAA (Lys) to TGG (Trp) has no stop-free ordering under vertebrate mito.

        The pair is still dropped from Dn/Ds (issue #90's chosen fix keeps that
        behavior), but GeneStats now counts how often it happens.
        """
        ingroup = write_vcf(tmp_path / "blocked_in", [])
        outgroup = write_vcf(
            tmp_path / "blocked_out",
            [
                "chr1\t7\t.\tA\tT\t30\tPASS\t.\tGT\t1/1",
                "chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
                "chr1\t9\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
            ],
            samples=OUTGROUP,
        )
        poly, stats = extract_gene_data(
            ingroup, outgroup, simple_cds, synthetic_ref, genetic_code=GeneticCode(table_id=2)
        )
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.stop_blocked_pairs == 1


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
    """A query that cannot run raises. It is not an empty result."""

    def test_unindexed_ingroup_raises(self, genome, ingroup_vcf_plain):
        with pytest.raises(VcfQueryError, match="index"):
            _extract(genome, "g_plus", vcf_path=ingroup_vcf_plain)

    def test_unindexed_outgroup_raises(self, genome, outgroup_vcf_plain):
        with pytest.raises(VcfQueryError, match="index"):
            _extract(genome, "g_plus", outgroup_vcf_path=outgroup_vcf_plain)


class TestProbeVcf:
    """probe_vcf proves a file answers region queries and reports its contigs."""

    def test_indexed_vcf_reports_its_contigs(self, genome):
        assert probe_vcf(genome.ingroup_vcf) == ["chr1", "chr2", "chrZ"]

    def test_plain_vcf_has_no_index_and_the_message_says_how_to_fix_it(self, ingroup_vcf_plain):
        with pytest.raises(VcfQueryError, match=r"has no usable index.*tabix"):
            probe_vcf(ingroup_vcf_plain)

    def test_bgzipped_vcf_without_its_index(self, genome, tmp_path):
        copy = tmp_path / "noindex.vcf.gz"
        copy.write_bytes(genome.ingroup_vcf.read_bytes())
        with pytest.raises(VcfQueryError, match="has no usable index"):
            probe_vcf(copy)

    def test_header_without_contig_lines_still_proves_the_index(self, tmp_path, write_vcf):
        """With no contig lines the header names nothing to query, and the proof must still run."""
        plain = write_vcf(tmp_path / "nocontig", [], contigs={}).with_suffix("")
        with pytest.raises(VcfQueryError, match="has no usable index"):
            probe_vcf(plain)

    def test_unreadable_file(self, not_a_vcf):
        with pytest.raises(VcfQueryError, match="cannot be opened"):
            probe_vcf(not_a_vcf)

    def test_htslib_index_error_goes_to_the_logger(self, ingroup_vcf_plain, caplog):
        """htslib's own complaint is captured at DEBUG instead of leaking to the terminal."""
        with caplog.at_level(logging.DEBUG, logger="mkado.io.vcf"):
            with pytest.raises(VcfQueryError):
                probe_vcf(ingroup_vcf_plain)
        assert any(r.message.startswith("htslib:") and "index" in r.message for r in caplog.records)


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


SYN_SITE = "chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t"
"""Codon 2, GCC to GCT (Ala to Ala). Append the sample genotypes, tab-separated."""


class TestCalledAlleleCounts:
    """Frequencies count each called allele and leave out each missing one."""

    def poly(self, tmp_path, write_vcf, synthetic_ref, simple_cds, genotypes, **kwargs):
        samples = [f"s{i}" for i in range(len(genotypes))]
        ingroup = write_vcf(tmp_path / "in", [SYN_SITE + "\t".join(genotypes)], samples=samples)
        poly, _ = extract_gene_data(ingroup, None, simple_cds, synthetic_ref, **kwargs)
        return poly.polymorphisms

    def test_partly_missing_ref_call_counts_one_ref(
        self, tmp_path, write_vcf, synthetic_ref, simple_cds
    ):
        """0/. adds one REF copy, not two."""
        polys = self.poly(
            tmp_path, write_vcf, synthetic_ref, simple_cds, ["0/0", "0/1", "0/.", "0/."]
        )
        _assert_polys(polys, [(1 / 6, "S")])

    def test_partly_missing_alt_call_counts_one_alt(
        self, tmp_path, write_vcf, synthetic_ref, simple_cds
    ):
        """./1 adds one ALT copy and no REF copy."""
        polys = self.poly(
            tmp_path, write_vcf, synthetic_ref, simple_cds, ["0/0", "0/1", "./1", "./1"]
        )
        _assert_polys(polys, [(3 / 6, "S")])

    def test_haploid_call_counts_one_copy(self, tmp_path, write_vcf, synthetic_ref, simple_cds):
        """A haploid call adds one copy, next to diploid calls in the same record."""
        polys = self.poly(tmp_path, write_vcf, synthetic_ref, simple_cds, ["1", "0/0", "0/0"])
        _assert_polys(polys, [(1 / 5, "S")])

    def test_haploid_singleton_is_removed(self, tmp_path, write_vcf, synthetic_ref, simple_cds):
        """One haploid ALT call is one derived copy, so --no-singletons removes it."""
        polys = self.poly(
            tmp_path,
            write_vcf,
            synthetic_ref,
            simple_cds,
            ["1", "0", "0", "0"],
            no_singletons=True,
        )
        assert polys == []

    def test_triploid_calls_count_every_allele(
        self, tmp_path, write_vcf, synthetic_ref, simple_cds
    ):
        """Each of the three alleles of a triploid call is counted."""
        polys = self.poly(tmp_path, write_vcf, synthetic_ref, simple_cds, ["0/0/1", "0/0/0"])
        _assert_polys(polys, [(1 / 6, "S")])

    def test_partly_missing_call_keeps_a_fixed_difference(
        self, tmp_path, write_vcf, synthetic_ref, simple_cds, outgroup_vcf_empty
    ):
        """./1 at a site fixed for ALT adds no REF copy, so the site stays fixed."""
        ingroup = write_vcf(
            tmp_path / "fixed_partial", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1\t1/1\t./1\t1/1"]
        )
        poly, _ = extract_gene_data(ingroup, outgroup_vcf_empty, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert (poly.dn, poly.ds) == (1, 0)


class TestOutgroupCalledAlleles:
    """The outgroup carries ALT only when every allele of its call is called ALT."""

    @pytest.mark.parametrize(
        "genotype, expected",
        [("1/.", (0, 0)), ("./1", (0, 0)), ("1", (1, 0)), ("1/1", (1, 0))],
    )
    def test_outgroup_call(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, genotype, expected
    ):
        """A partly missing call is unresolved and reads as the reference base."""
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", [f"chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t{genotype}"], samples=OUTGROUP
        )
        poly, _ = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == expected


@pytest.fixture(params=["N", "R"])
def ambiguous_ref(request, tmp_path: Path) -> Path:
    """The synthetic reference with an ambiguous base in codon 3: ATG GCC A?A."""
    ref_path = tmp_path / "ambiguous_ref.fa"
    ref_path.write_text(f">chr1\nATGGCCA{request.param}ATTTTTTTTTTTTTTTTTTTTTT\n")
    pysam.faidx(str(ref_path))
    return ref_path


class TestAmbiguousReference:
    """A codon with an ambiguous reference base is skipped and counted as ambiguous."""

    def test_divergence_is_not_reported_as_stop_blocked(
        self, ambiguous_ref, simple_cds, write_vcf, tmp_path
    ):
        """An outgroup ALT in the ambiguous codon is skipped, not counted as stop-blocked."""
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t9\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, ambiguous_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.stop_blocked_pairs == 0
        assert stats.skipped_ambiguous == 1

    def test_polymorphism_is_counted_as_skipped(
        self, ambiguous_ref, simple_cds, write_vcf, tmp_path
    ):
        """An ingroup SNP in the ambiguous codon is skipped and counted."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t9\t.\tA\tG\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0"]
        )
        poly, stats = extract_gene_data(ingroup, None, simple_cds, ambiguous_ref)
        assert poly.polymorphisms == []
        assert stats.skipped_ambiguous == 1


class TestSpanningDeletion:
    """A ``*`` ALT is a deletion that spans the site, so it is skipped as an indel."""

    def test_ingroup_star_is_counted_as_indel(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t6\t.\tC\t*\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0"]
        )
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert stats.skipped_indels == 1
        assert stats.skipped_ambiguous == 0

    def test_outgroup_star_does_not_void_a_difference(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """The deletion reads as the reference base, so the A>G difference in codon 3 stays."""
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out",
            [
                "chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
                "chr1\t9\t.\tA\t*\t30\tPASS\t.\tGT\t1/1",
            ],
            samples=OUTGROUP,
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (1, 0)
        assert stats.skipped_ambiguous == 0


class TestCdsPastContigEnd:
    def test_cds_past_contig_end_is_an_error(self, synthetic_ref, write_vcf, tmp_path):
        """A CDS that ends past its contig does not match the reference, so the gene fails."""
        cds = CdsRegion(
            gene_id="long", transcript_id="tx", chrom="chr1", exons=[(0, 33)], strand="+"
        )
        ingroup = write_vcf(tmp_path / "in", [])
        with pytest.raises(ValueError, match="past the end of chr1"):
            extract_gene_data(ingroup, None, cds, synthetic_ref)


class TestRecordWithoutCalledAlt:
    """A record whose called ingroup alleles are all REF is the same as no record."""

    def extract(self, synthetic_ref, simple_cds, write_vcf, tmp_path, genotypes):
        """Ingroup record in codon 3 at position 7; outgroup fixed for G at position 8."""
        ingroup = write_vcf(tmp_path / "in", [f"chr1\t7\t.\tA\tC\t30\tPASS\t.\tGT\t{genotypes}"])
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        return extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)

    @pytest.mark.parametrize("genotypes", ["0/0\t0/0\t0/0\t0/0", "0/0\t./.\t0/0\t./."])
    def test_keeps_the_fixed_difference_in_its_codon(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, genotypes
    ):
        poly, stats = self.extract(synthetic_ref, simple_cds, write_vcf, tmp_path, genotypes)
        assert (poly.dn, poly.ds) == (1, 0)
        assert poly.polymorphisms == []
        assert stats == GeneStats()

    def test_real_polymorphism_still_voids_the_codon(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """A called ALT allele in the same codon is a polymorphism, so no fixed difference."""
        poly, _ = self.extract(synthetic_ref, simple_cds, write_vcf, tmp_path, "0/0\t0/1\t0/0\t0/0")
        assert (poly.dn, poly.ds) == (0, 0)
        assert len(poly.polymorphisms) == 1


SITE_6 = "chr1\t6\t.\tC\tT\t30\t{filter}\t.\tGT\t0/0\t0/1\t0/0\t0/0"
"""A synonymous SNP in codon 2 (GCC to GCT) with the FILTER value left open."""


class TestInvariantRecords:
    def test_dot_alt_is_not_counted(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        """A record with no ALT is an invariant site, not a multi-allelic one."""
        ingroup = write_vcf(
            tmp_path / "in",
            [
                "chr1\t5\t.\tC\t.\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0",
                "chr1\t6\t.\tC\t.\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0",
            ],
        )
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert stats == GeneStats()


class TestFilterColumn:
    def test_filtered_record_is_skipped_and_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        ingroup = write_vcf(tmp_path / "in", [SITE_6.format(filter="LowQual")])
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert stats.skipped_filtered == 1

    def test_include_filtered_uses_the_record(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        ingroup = write_vcf(tmp_path / "in", [SITE_6.format(filter="LowQual")])
        poly, stats = extract_gene_data(
            ingroup, None, simple_cds, synthetic_ref, include_filtered=True
        )
        _assert_polys(poly.polymorphisms, [(0.125, "S")])
        assert stats.skipped_filtered == 0

    @pytest.mark.parametrize("value", ["PASS", "."])
    def test_pass_and_missing_filter_are_used(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, value
    ):
        ingroup = write_vcf(tmp_path / "in", [SITE_6.format(filter=value)])
        poly, _ = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        _assert_polys(poly.polymorphisms, [(0.125, "S")])

    def test_filtered_outgroup_record_is_skipped_and_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tA\tG\t30\tLowQual\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.outgroup_skipped == 1


class TestReferenceMismatch:
    def test_ingroup_ref_mismatch_is_skipped_and_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """The VCF says REF T at position 6, but the FASTA base there is C."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t6\t.\tT\tA\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0"]
        )
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert stats.skipped_ref_mismatch == 1

    def test_outgroup_ref_mismatch_is_counted(self, synthetic_ref, simple_cds, write_vcf, tmp_path):
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t6\t.\tT\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.outgroup_skipped == 1


class TestOutgroupSkipCount:
    @pytest.mark.parametrize(
        "record",
        [
            "chr1\t5\t.\tCC\tC\t30\tPASS\t.\tGT\t1/1",
            "chr1\t8\t.\tA\tG,T\t30\tPASS\t.\tGT\t1/1",
        ],
        ids=["indel", "multi_allelic"],
    )
    def test_skipped_outgroup_record_is_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, record
    ):
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(tmp_path / "out", [record], samples=OUTGROUP)
        _, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert stats.outgroup_skipped == 1

    def test_outgroup_record_without_genotype_is_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(tmp_path / "out", ["chr1\t8\t.\tA\tG\t30\tPASS\t."], samples=[])
        _, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert stats.outgroup_skipped == 1


class TestGenotypedSamples:
    def test_names_of_a_genotyped_vcf(self, write_vcf, tmp_path):
        assert genotyped_samples(write_vcf(tmp_path / "in", [])) == ["s1", "s2", "s3", "s4"]

    def test_sites_only_vcf_has_none(self, write_vcf, tmp_path):
        assert genotyped_samples(write_vcf(tmp_path / "sites", [], samples=[])) == []

    def test_vcf_without_gt_field_has_none(self, tmp_path):
        """Sample columns that hold only DP give no genotypes."""
        plain = tmp_path / "dp.vcf"
        plain.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=31>\n"
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
            "chr1\t6\t.\tC\tT\t30\tPASS\t.\tDP\t5\n"
        )
        gz = tmp_path / "dp.vcf.gz"
        pysam.tabix_compress(str(plain), str(gz), force=True)
        pysam.tabix_index(str(gz), preset="vcf", force=True)
        assert genotyped_samples(gz) == []


class TestReferenceMismatchCodon:
    """The FASTA base at a REF mismatch cannot be trusted, so its codon is not used."""

    def test_ingroup_mismatch_removes_the_codon_from_divergence(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """REF T where the FASTA has A at position 7, in codon 3 with an outgroup difference."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t7\t.\tT\tC\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0"]
        )
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.skipped_ref_mismatch == 1

    def test_snp_sharing_a_codon_with_a_mismatch_is_skipped(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """The C>T SNP at position 6 shares codon 2 with the mismatched record at position 5."""
        ingroup = write_vcf(
            tmp_path / "in",
            [
                "chr1\t5\t.\tT\tG\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0",
                "chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0",
            ],
        )
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert stats.skipped_ref_mismatch == 2

    def test_outgroup_mismatch_removes_the_codon_from_divergence(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """REF G where the FASTA has A at position 9, in codon 3 with the A>G difference."""
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out",
            [
                "chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1",
                "chr1\t9\t.\tG\tC\t30\tPASS\t.\tGT\t1/1",
            ],
            samples=OUTGROUP,
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        # The mismatched record and the A>G difference that shares its codon.
        assert stats.outgroup_skipped == 2

    def test_filtered_mismatch_still_removes_its_codon(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """A filtered record with a REF mismatch still shows the FASTA base is wrong."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t7\t.\tT\tC\t30\tLowQual\t.\tGT\t0/0\t0/1\t0/0\t0/0"]
        )
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.skipped_ref_mismatch == 1
        assert stats.skipped_filtered == 0

    def test_fixed_alt_snp_in_a_mismatch_codon_is_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, outgroup_vcf_empty
    ):
        """A site fixed for ALT in a codon with a REF mismatch is dropped and counted."""
        ingroup = write_vcf(
            tmp_path / "in",
            [
                "chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1",
                "chr1\t9\t.\tG\tC\t30\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0",
            ],
        )
        poly, stats = extract_gene_data(ingroup, outgroup_vcf_empty, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.skipped_ref_mismatch == 2


class TestSwappedRefAlt:
    """A VCF REF that is not the FASTA base, with an ALT that is, is read with the two swapped."""

    def test_ingroup_swap_gives_a_polymorphism(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """REF T, ALT C at position 6 (FASTA C): three copies of T in eight."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t6\t.\tT\tC\t30\tPASS\t.\tGT\t0/0\t0/1\t1/1\t1/1"]
        )
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        _assert_polys(poly.polymorphisms, [(0.375, "S")])
        assert stats.swapped_ref_alt == 1
        assert stats.skipped_ref_mismatch == 0

    def test_ingroup_swap_can_be_fixed(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, outgroup_vcf_empty
    ):
        """REF G, ALT A at position 8 (FASTA A), all 0/0: the ingroup is fixed for G."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t8\t.\tG\tA\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0"]
        )
        poly, _ = extract_gene_data(ingroup, outgroup_vcf_empty, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (1, 0)

    def test_outgroup_swap_with_the_non_reference_call(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """A 0/0 call on a swapped record means the outgroup carries G, not the FASTA A."""
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tG\tA\t30\tPASS\t.\tGT\t0/0"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (1, 0)
        assert stats.swapped_ref_alt == 1

    def test_outgroup_swap_with_the_reference_call(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """A 1/1 call on a swapped record means the outgroup carries the FASTA base."""
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tG\tA\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.outgroup_skipped == 0

    def test_filtered_swap_does_not_remove_its_codon(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path
    ):
        """A swapped record agrees with the FASTA, so only the FILTER rule skips it."""
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t9\t.\tG\tA\t30\tLowQual\t.\tGT\t0/0\t0/1\t0/0\t0/0"]
        )
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (1, 0)
        assert stats.skipped_filtered == 1
        assert stats.skipped_ref_mismatch == 0


class TestUnresolvedOutgroupCall:
    """An outgroup call that cannot resolve one base reads as REF, and is counted."""

    @pytest.mark.parametrize("genotype", ["./.", "1/.", "./1", "0/1"])
    def test_unresolved_call_is_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, genotype
    ):
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", [f"chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t{genotype}"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.outgroup_skipped == 1

    @pytest.mark.parametrize("genotype, dn", [("0/0", 0), ("1/1", 1)])
    def test_resolved_call_is_not_counted(
        self, synthetic_ref, simple_cds, write_vcf, tmp_path, genotype, dn
    ):
        ingroup = write_vcf(tmp_path / "in", [])
        outgroup = write_vcf(
            tmp_path / "out", [f"chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t{genotype}"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert poly.dn == dn
        assert stats.outgroup_skipped == 0


# Position 7 (1-based) is the first base of AAA, so A to T makes the stop codon TAA.
AAA_TO_TAA = "chr1\t7\t.\tA\tT\t30\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0"


class TestStopCodons:
    """A change to or from a stop codon is skipped and counted, in both polymorphism and divergence."""

    def test_snp_that_makes_a_stop_is_skipped(self, tmp_path, write_vcf, synthetic_ref, simple_cds):
        ingroup = write_vcf(tmp_path / "in_stop", [AAA_TO_TAA])
        poly, stats = extract_gene_data(ingroup, None, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert stats.skipped_stop == 1

    def test_fixed_difference_to_a_stop_is_counted(
        self, tmp_path, write_vcf, synthetic_ref, simple_cds
    ):
        ingroup = write_vcf(tmp_path / "in_empty", [])
        outgroup = write_vcf(
            tmp_path / "out_stop",
            ["chr1\t7\t.\tA\tT\t30\tPASS\t.\tGT\t1/1"],
            samples=OUTGROUP,
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert (poly.dn, poly.ds) == (0, 0)
        assert stats.skipped_stop == 1

    def test_snp_that_makes_a_stop_does_not_void_the_codon(
        self, tmp_path, write_vcf, synthetic_ref, simple_cds
    ):
        """The stop allele is missing data, so the ingroup codon is AAA, as in the FASTA path."""
        ingroup = write_vcf(tmp_path / "in_stop", [AAA_TO_TAA])
        # The outgroup is fixed for AGA (Arg) in the same codon: one replacement from AAA (Lys).
        outgroup = write_vcf(
            tmp_path / "out_arg", ["chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, simple_cds, synthetic_ref)
        assert poly.polymorphisms == []
        assert (poly.dn, poly.ds) == (1, 0)
        assert stats.skipped_stop == 1

    def test_ref_stop_with_non_stop_alt_matches_fasta(self, tmp_path, write_vcf):
        """The REF stop carriers are missing data, so the ingroup codon is the ALT codon.

        REF codon 1 is TAA. One ingroup copy carries TCA (Ser); the outgroup is fixed for
        CAA (Gln). TCA to CAA is two replacements, the same as the FASTA path counts.
        """
        from mkado.analysis.mk_test import mk_test
        from tests.builders import sequence_set

        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\nATGTAAGCCTTTTTTTTTTTTTTTTTTTTTTT\n")
        pysam.faidx(str(ref))
        cds = CdsRegion(gene_id="g", transcript_id="t", chrom="chr1", exons=[(0, 9)], strand="+")
        ingroup = write_vcf(
            tmp_path / "in", ["chr1\t5\t.\tA\tC\t30\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0"]
        )
        outgroup = write_vcf(
            tmp_path / "out", ["chr1\t4\t.\tT\tC\t30\tPASS\t.\tGT\t1/1"], samples=OUTGROUP
        )
        poly, stats = extract_gene_data(ingroup, outgroup, cds, ref)

        fasta = mk_test(
            sequence_set(["ATGTAAGCC"] * 7 + ["ATGTCAGCC"]), sequence_set(["ATGCAAGCC"])
        )
        assert (poly.dn, poly.ds) == (fasta.dn, fasta.ds) == (2, 0)
        assert poly.polymorphisms == []
        assert stats.skipped_stop == 1

    def test_warning_counts_stop_skips(self):
        from mkado.vcf_workers import _skip_warning

        assert "2 with stop codons" in _skip_warning("g", GeneStats(skipped_stop=2))
