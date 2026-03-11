"""Tests for VCF input functionality.

These tests require cyvcf2 and pysam. They are skipped if not installed.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

pytest.importorskip("cyvcf2", reason="cyvcf2 not installed")
pytest.importorskip("pysam", reason="pysam not installed")

import cyvcf2  # noqa: E402
import pysam  # noqa: E402

from mkado.core.cds import CdsRegion  # noqa: E402
from mkado.io.gff import parse_gff3  # noqa: E402
from mkado.io.vcf import (  # noqa: E402
    _check_vcf_deps,
    _reconstruct_codon_with_sub,
    extract_gene_data,
)


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


def _write_vcf(path: Path, header_lines: list[str], records: list[str]) -> Path:
    """Write a VCF file and create tabix index."""
    import subprocess

    vcf_path = path.with_suffix(".vcf")
    lines = []
    lines.append("##fileformat=VCFv4.2")
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append(
        "##contig=<ID=chr1,length=30>"
    )
    for h in header_lines:
        lines.append(h)
    lines.extend(records)

    vcf_path.write_text("\n".join(lines) + "\n")

    # bgzip and tabix
    bgz_path = path.with_suffix(".vcf.gz")
    subprocess.run(["bgzip", "-c", str(vcf_path)], stdout=open(str(bgz_path), "wb"), check=True)
    subprocess.run(["tabix", "-p", "vcf", str(bgz_path)], check=True)

    return bgz_path


@pytest.fixture
def ingroup_vcf_synonymous(tmp_path: Path) -> Path:
    """Ingroup VCF with one synonymous polymorphism.

    Site at position 5 (0-based): C->T in codon GCC -> GCT (both Ala).
    4 diploid samples: 3 hom-ref, 1 het => alt_freq = 1/8 = 0.125
    VCF is 1-based, so pos=6.
    """
    header = [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsamp1\tsamp2\tsamp3\tsamp4"
    ]
    records = [
        "chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1"
    ]
    return _write_vcf(tmp_path / "ingroup_syn", header, records)


@pytest.fixture
def ingroup_vcf_nonsyn(tmp_path: Path) -> Path:
    """Ingroup VCF with one nonsynonymous polymorphism.

    Site at position 3 (0-based): G->A in codon GCC -> ACC (Ala -> Thr).
    4 diploid samples: 2 hom-ref, 2 het => alt_freq = 2/8 = 0.25
    VCF is 1-based, so pos=4.
    """
    header = [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsamp1\tsamp2\tsamp3\tsamp4"
    ]
    records = [
        "chr1\t4\t.\tG\tA\t30\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1"
    ]
    return _write_vcf(tmp_path / "ingroup_nonsyn", header, records)


@pytest.fixture
def outgroup_vcf_divergent(tmp_path: Path) -> Path:
    """Outgroup VCF with a fixed difference at codon 3 (AAA -> AGA = Lys -> Arg).

    Site at position 7 (0-based): A->G. VCF pos=8.
    Single sample, hom-alt.
    """
    header = [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\toutgroup"
    ]
    records = [
        "chr1\t8\t.\tA\tG\t30\tPASS\t.\tGT\t1/1"
    ]
    return _write_vcf(tmp_path / "outgroup", header, records)


@pytest.fixture
def outgroup_vcf_empty(tmp_path: Path) -> Path:
    """Outgroup VCF with no variants (all same as reference)."""
    header = [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\toutgroup"
    ]
    return _write_vcf(tmp_path / "outgroup_empty", header, [])


# ---- Tests ----


class TestCheckDeps:
    def test_deps_available(self):
        """cyvcf2 and pysam should be available in test environment."""
        _check_vcf_deps()  # Should not raise


class TestCodonReconstruction:
    def test_plus_strand_substitution(self, synthetic_ref, simple_cds):
        fasta = pysam.FastaFile(str(synthetic_ref))
        ref_fetch = lambda chrom, pos: fasta.fetch(chrom, pos, pos + 1).upper()

        # Substitute position 5 (C->T) in codon GCC (index 1)
        ref_codon, alt_codon = _reconstruct_codon_with_sub(
            simple_cds, 1, ref_fetch, 5, "T"
        )
        assert ref_codon == "GCC"
        assert alt_codon == "GCT"
        fasta.close()

    def test_nonsynonymous_substitution(self, synthetic_ref, simple_cds):
        fasta = pysam.FastaFile(str(synthetic_ref))
        ref_fetch = lambda chrom, pos: fasta.fetch(chrom, pos, pos + 1).upper()

        # Substitute position 3 (G->A) in codon GCC (index 1)
        ref_codon, alt_codon = _reconstruct_codon_with_sub(
            simple_cds, 1, ref_fetch, 3, "A"
        )
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

    def test_gene_id_propagated(self, synthetic_ref, simple_cds, ingroup_vcf_synonymous, outgroup_vcf_empty):
        poly_data, _ = extract_gene_data(
            vcf_path=ingroup_vcf_synonymous,
            outgroup_vcf_path=outgroup_vcf_empty,
            cds=simple_cds,
            ref_fasta_path=synthetic_ref,
        )
        assert poly_data.gene_id == "test_gene"


class TestPolarization:
    def test_polarization_flips_frequency(
        self, synthetic_ref, simple_cds, tmp_path
    ):
        """If outgroup carries ALT, derived freq should be 1 - alt_freq."""
        # Create ingroup with a SNP at pos 5 (C->T, synonymous)
        header = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4"
        ]
        ingroup_records = [
            "chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1"
        ]
        ingroup_vcf = _write_vcf(tmp_path / "ig_polar", header, ingroup_records)

        # Outgroup also carries T at this position
        out_header = [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\toutgroup"
        ]
        out_records = [
            "chr1\t6\t.\tC\tT\t30\tPASS\t.\tGT\t1/1"
        ]
        outgroup_vcf = _write_vcf(tmp_path / "og_polar", out_header, out_records)

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
