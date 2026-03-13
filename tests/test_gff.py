"""Tests for GFF3 parser."""

from __future__ import annotations

import gzip
import textwrap
from pathlib import Path

import pytest

from mkado.io.gff import parse_gff3


@pytest.fixture
def simple_gff3(tmp_path: Path) -> Path:
    """Create a simple GFF3 file with two genes (CDS lengths divisible by 3).

    GFF3 coords are 1-based inclusive, so start=101 end=109 = 9 bases.
    Converted to 0-based half-open: (100, 109).
    """
    # GeneA: CDS exon1 = 101-109 (9bp) + exon2 = 301-309 (9bp) = 18bp = 6 codons
    # GeneB: CDS = 501-509 (9bp) = 3 codons
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t101\t400\t.\t+\t.\tID=gene1;Name=GeneA
        chr1\t.\tmRNA\t101\t400\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t101\t109\t.\t+\t0\tID=cds1;Parent=tx1
        chr1\t.\tCDS\t301\t309\t.\t+\t0\tID=cds2;Parent=tx1
        chr2\t.\tgene\t501\t600\t.\t-\t.\tID=gene2;Name=GeneB
        chr2\t.\tmRNA\t501\t600\t.\t-\t.\tID=tx2;Parent=gene2
        chr2\t.\tCDS\t501\t509\t.\t-\t0\tID=cds3;Parent=tx2
    """)
    gff_path = tmp_path / "test.gff3"
    gff_path.write_text(content)
    return gff_path


@pytest.fixture
def multi_transcript_gff3(tmp_path: Path) -> Path:
    """GFF3 with a gene that has two transcripts (longest should be selected)."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t101\t600\t.\t+\t.\tID=gene1;Name=GeneA
        chr1\t.\tmRNA\t101\t300\t.\t+\t.\tID=tx_short;Parent=gene1
        chr1\t.\tCDS\t101\t199\t.\t+\t0\tID=cds_s1;Parent=tx_short
        chr1\t.\tmRNA\t101\t600\t.\t+\t.\tID=tx_long;Parent=gene1
        chr1\t.\tCDS\t101\t199\t.\t+\t0\tID=cds_l1;Parent=tx_long
        chr1\t.\tCDS\t301\t501\t.\t+\t0\tID=cds_l2;Parent=tx_long
    """)
    gff_path = tmp_path / "multi_tx.gff3"
    gff_path.write_text(content)
    return gff_path


def test_parse_simple_gff3(simple_gff3: Path):
    regions = parse_gff3(simple_gff3)
    assert len(regions) == 2

    # Find gene1 (plus strand, two exons)
    gene_a = [r for r in regions if r.gene_id == "gene1"]
    assert len(gene_a) == 1
    cds = gene_a[0]
    assert cds.chrom == "chr1"
    assert cds.strand == "+"
    assert len(cds.exons) == 2
    # GFF3 coordinates are 1-based inclusive, converted to 0-based half-open
    assert cds.exons[0] == (100, 109)
    assert cds.exons[1] == (300, 309)
    assert cds.cds_length == 18  # 9 + 9
    assert cds.is_valid()

    # gene2 (minus strand, single exon)
    gene_b = [r for r in regions if r.gene_id == "gene2"]
    assert len(gene_b) == 1
    cds_b = gene_b[0]
    assert cds_b.chrom == "chr2"
    assert cds_b.strand == "-"
    assert cds_b.cds_length == 9


def test_parse_filters_invalid_cds(simple_gff3: Path):
    """parse_gff3 should skip genes where CDS length % 3 != 0."""
    regions = parse_gff3(simple_gff3)
    # Both genes have CDS lengths not divisible by 3 (200, 100)
    # They should be filtered out
    assert all(r.is_valid() for r in regions)


def test_parse_with_gene_filter(simple_gff3: Path):
    regions = parse_gff3(simple_gff3, gene_ids={"GeneA"})
    gene_ids = {r.gene_id for r in regions}
    assert "gene2" not in gene_ids


def test_selects_longest_transcript(multi_transcript_gff3: Path):
    regions = parse_gff3(multi_transcript_gff3)
    if regions:
        gene_a = [r for r in regions if r.gene_id == "gene1"]
        if gene_a:
            cds = gene_a[0]
            # tx_long has 99 + 201 = 300 bases (longer than tx_short's 99)
            assert cds.transcript_id == "tx_long"


@pytest.fixture
def valid_gff3(tmp_path: Path) -> Path:
    """GFF3 where CDS lengths are divisible by 3."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene1;Name=ValidGene
        chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "valid.gff3"
    gff_path.write_text(content)
    return gff_path


def test_parse_valid_cds(valid_gff3: Path):
    regions = parse_gff3(valid_gff3)
    assert len(regions) == 1
    cds = regions[0]
    assert cds.gene_id == "gene1"
    assert cds.cds_length == 9
    assert cds.num_codons() == 3
    assert cds.is_valid()


def test_empty_file(tmp_path: Path):
    gff_path = tmp_path / "empty.gff3"
    gff_path.write_text("##gff-version 3\n")
    regions = parse_gff3(gff_path)
    assert regions == []


def test_comments_and_blank_lines(tmp_path: Path):
    content = textwrap.dedent("""\
        ##gff-version 3
        # This is a comment

        chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene1;Name=G1
        chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "comments.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1


def test_parse_gzipped_gff3(tmp_path: Path):
    """GFF3 files compressed with gzip should be parsed correctly."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene1;Name=GzGene
        chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "test.gff3.gz"
    with gzip.open(gff_path, "wt") as f:
        f.write(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"
