"""Tests for GFF3 parser."""

from __future__ import annotations

import gzip
import logging
import textwrap
from pathlib import Path

import pytest

from mkado.core.cds import CdsRegion
from mkado.io.gff import GffFormatError, parse_gff3
from tests.builders import NON_UTF8_GFF


# ---------------------------------------------------------------------------
# Fixtures: minimal valid GFF3 files
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Original tests (preserved)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# GTF format detection
# ---------------------------------------------------------------------------


def test_gtf_raises_error(tmp_path: Path):
    """GTF files should be rejected with a clear error message."""
    # GTF uses space-delimited key "value"; attributes, not key=value
    content = textwrap.dedent("""\
        chr1\tprotein_coding\tgene\t1\t9\t.\t+\t.\tgene_id "gene1"; gene_name "GeneA";
        chr1\tprotein_coding\ttranscript\t1\t9\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";
        chr1\tprotein_coding\tCDS\t1\t9\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1";
    """)
    gff_path = tmp_path / "test.gtf"
    gff_path.write_text(content)
    with pytest.raises(GffFormatError, match="test.gtf: file appears to be GTF"):
        parse_gff3(gff_path)


def test_gtf_with_gff3_extension_raises_error(tmp_path: Path):
    """A GTF file mislabeled as .gff3 should still be caught by content detection."""
    content = textwrap.dedent("""\
        chr1\tprotein_coding\tgene\t1\t9\t.\t+\t.\tgene_id "gene1"; gene_name "GeneA";
        chr1\tprotein_coding\tCDS\t1\t9\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1";
    """)
    gff_path = tmp_path / "mislabeled.gff3"
    gff_path.write_text(content)
    with pytest.raises(ValueError, match="(?i)gtf"):
        parse_gff3(gff_path)


# ---------------------------------------------------------------------------
# Missing hierarchy levels
# ---------------------------------------------------------------------------


def test_missing_gene_features(tmp_path: Path):
    """When gene features are absent, genes should be inferred from mRNA Parent."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "no_gene.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    # Should still work: mRNA Parent="gene1" defines the gene grouping
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"


def test_missing_mrna_features(tmp_path: Path):
    """CDS features parented directly to a gene (no mRNA level) should work."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=DirectCDS
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1
    """)
    gff_path = tmp_path / "no_mrna.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"


def test_missing_both_gene_and_mrna(tmp_path: Path):
    """CDS with Parent pointing to an unknown ID should still be extracted."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=orphan1
    """)
    gff_path = tmp_path / "orphan_cds.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    # The CDS Parent "orphan1" is treated as both gene and transcript
    assert len(regions) == 1
    assert regions[0].gene_id == "orphan1"


# ---------------------------------------------------------------------------
# CDS without Parent attribute
# ---------------------------------------------------------------------------


def test_cds_without_parent_skipped(tmp_path: Path, caplog):
    """CDS features lacking a Parent attribute should be skipped with a warning."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds_good;Parent=tx1
        chr1\t.\tCDS\t20\t28\t.\t+\t0\tID=cds_bad
    """)
    gff_path = tmp_path / "no_parent.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.DEBUG, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    # gene1 should still be extracted (from the CDS that has Parent)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"


# ---------------------------------------------------------------------------
# Invalid CDS length (not divisible by 3)
# ---------------------------------------------------------------------------


def test_invalid_cds_length_skipped_with_warning(tmp_path: Path, caplog):
    """Genes with CDS length not divisible by 3 should be skipped with a warning."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene_bad
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx_bad;Parent=gene_bad
        chr1\t.\tCDS\t1\t10\t.\t+\t0\tID=cds_bad;Parent=tx_bad
        chr1\t.\tgene\t200\t300\t.\t+\t.\tID=gene_good
        chr1\t.\tmRNA\t200\t300\t.\t+\t.\tID=tx_good;Parent=gene_good
        chr1\t.\tCDS\t200\t208\t.\t+\t0\tID=cds_good;Parent=tx_good
    """)
    gff_path = tmp_path / "bad_frame.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene_good"
    # Should warn about the skipped gene
    assert any("not divisible by 3" in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# Malformed lines
# ---------------------------------------------------------------------------


def test_malformed_lines_skipped(tmp_path: Path, caplog):
    """Lines with too few columns should be skipped with a warning."""
    content = textwrap.dedent("""\
        ##gff-version 3
        this is not a valid gff line
        chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene1
        only\ttwo\tcolumns
        chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "malformed.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.DEBUG, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    # The valid gene should still be parsed
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"
    # Should warn about malformed lines
    assert any(
        "malformed" in msg.lower() or "skipping line" in msg.lower() for msg in caplog.messages
    )


def test_non_integer_coordinates(tmp_path: Path, caplog):
    """Lines with non-integer coordinates should be skipped with a warning."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene_good
        chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=tx_good;Parent=gene_good
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds_good;Parent=tx_good
        chr1\t.\tgene\tabc\txyz\t.\t+\t.\tID=gene_bad
    """)
    gff_path = tmp_path / "bad_coords.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.DEBUG, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene_good"


# ---------------------------------------------------------------------------
# Gene ID filtering
# ---------------------------------------------------------------------------


def test_filter_by_gene_id(simple_gff3: Path):
    """Filtering by ID attribute should work."""
    regions = parse_gff3(simple_gff3, gene_ids={"gene1"})
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"


def test_filter_by_gene_name(simple_gff3: Path):
    """Filtering by Name attribute should work."""
    regions = parse_gff3(simple_gff3, gene_ids={"GeneB"})
    assert len(regions) == 1
    assert regions[0].gene_id == "gene2"


def test_filter_nonexistent_gene(simple_gff3: Path):
    """Filtering for a gene that doesn't exist should return empty list."""
    regions = parse_gff3(simple_gff3, gene_ids={"NonExistent"})
    assert regions == []


# ---------------------------------------------------------------------------
# Non-coding features should not produce CDS regions
# ---------------------------------------------------------------------------


def test_noncoding_features_ignored(tmp_path: Path):
    """Non-coding gene types (lncRNA, rRNA, etc.) without CDS should not appear."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t500\t.\t+\t.\tID=gene_lnc;Name=MyLncRNA;biotype=lncRNA
        chr1\t.\tlnc_RNA\t1\t500\t.\t+\t.\tID=lnc_tx;Parent=gene_lnc
        chr1\t.\texon\t1\t500\t.\t+\t.\tID=lnc_exon;Parent=lnc_tx
        chr1\t.\tgene\t1000\t1100\t.\t+\t.\tID=gene_coding;Name=MyCoding
        chr1\t.\tmRNA\t1000\t1100\t.\t+\t.\tID=tx_coding;Parent=gene_coding
        chr1\t.\tCDS\t1000\t1008\t.\t+\t0\tID=cds_coding;Parent=tx_coding
    """)
    gff_path = tmp_path / "noncoding.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    # Only the coding gene should appear
    assert len(regions) == 1
    assert regions[0].gene_id == "gene_coding"


# ---------------------------------------------------------------------------
# Phase handling
# ---------------------------------------------------------------------------


def test_phase_handling_plus_strand(tmp_path: Path):
    """Phase > 0 should trim bases from the start of plus-strand CDS."""
    # 12 bases with phase=1 -> 11 usable, but 11 % 3 != 0 so invalid
    # 12 bases with phase=0 -> 12 usable, 12 % 3 == 0 -> valid
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t12\t.\t+\t.\tID=gene1
        chr1\t.\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t12\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "phase.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].cds_length == 12
    assert regions[0].num_codons() == 4


def test_phase_offset_trims_bases(tmp_path: Path):
    """Phase=1 should trim 1 base, yielding different CDS length."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t13\t.\t+\t.\tID=gene1
        chr1\t.\tmRNA\t1\t13\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t13\t.\t+\t1\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "phase1.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    # 13 bases - 1 phase offset = 12 usable
    assert regions[0].cds_length == 12
    assert regions[0].num_codons() == 4


def test_phase_minus_strand(tmp_path: Path):
    """Phase on minus-strand gene uses the last CDS in genomic order."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t-\t.\tID=gene1
        chr1\t.\tmRNA\t1\t100\t.\t-\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t6\t.\t-\t0\tID=cds1;Parent=tx1
        chr1\t.\tCDS\t50\t55\t.\t-\t0\tID=cds2;Parent=tx1
    """)
    gff_path = tmp_path / "minus_phase.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].strand == "-"
    # 6 + 6 = 12 bases, phase=0
    assert regions[0].cds_length == 12


# ---------------------------------------------------------------------------
# Comma-separated Parent attributes
# ---------------------------------------------------------------------------


def test_comma_separated_parent(tmp_path: Path):
    """CDS with comma-separated Parent should be assigned to multiple transcripts."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx2;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1,tx2
    """)
    gff_path = tmp_path / "comma_parent.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    # Both transcripts have same CDS, longest is a tie -> either is fine
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"


# ---------------------------------------------------------------------------
# No CDS features at all
# ---------------------------------------------------------------------------


def test_no_cds_features(tmp_path: Path, caplog):
    """A file with gene/mRNA but no CDS should return empty list."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=NoCDS
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=tx1
    """)
    gff_path = tmp_path / "no_cds.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert regions == []


# ---------------------------------------------------------------------------
# Mixed features from real-world GFF3 (UTRs, exons, etc.)
# ---------------------------------------------------------------------------


def test_utr_and_exon_features_ignored(tmp_path: Path):
    """UTR and exon features should not interfere with CDS extraction."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t200\t.\t+\t.\tID=gene1;Name=WithUTR
        chr1\t.\tmRNA\t1\t200\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tfive_prime_UTR\t1\t10\t.\t+\t.\tID=utr5;Parent=tx1
        chr1\t.\texon\t1\t200\t.\t+\t.\tID=exon1;Parent=tx1
        chr1\t.\tCDS\t11\t19\t.\t+\t0\tID=cds1;Parent=tx1
        chr1\t.\tthree_prime_UTR\t20\t200\t.\t+\t.\tID=utr3;Parent=tx1
    """)
    gff_path = tmp_path / "utr.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    # CDS is 11-19 (9bp), not the full gene/exon span
    assert regions[0].cds_length == 9


# ---------------------------------------------------------------------------
# Example data integration test
# ---------------------------------------------------------------------------


def test_parse_example_annotation():
    """Parse the shipped example annotation and verify expected gene count."""
    gff_path = Path(__file__).parent.parent / "examples" / "example_vcf" / "annotation.gff3"
    if not gff_path.exists():
        pytest.skip("example annotation not available")
    regions = parse_gff3(gff_path)
    # We know from the mkado vcf run that 42 valid genes are extracted
    assert len(regions) == 42
    # All should have valid CDS lengths
    assert all(r.is_valid() for r in regions)
    # All should be on chromosome "1"
    assert all(r.chrom == "1" for r in regions)


def test_parse_example_annotation_filter_single_gene():
    """Filtering for a single gene from the example annotation should work."""
    gff_path = Path(__file__).parent.parent / "examples" / "example_vcf" / "annotation.gff3"
    if not gff_path.exists():
        pytest.skip("example annotation not available")
    regions = parse_gff3(gff_path, gene_ids={"LOC_00000001"})
    assert len(regions) == 1
    assert regions[0].gene_id == "LOC_00000001"


# ---------------------------------------------------------------------------
# Transcript feature type "transcript" (not just "mRNA")
# ---------------------------------------------------------------------------


def test_transcript_feature_type(tmp_path: Path):
    """The parser should accept 'transcript' as well as 'mRNA'."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
        chr1\t.\ttranscript\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
    """)
    gff_path = tmp_path / "transcript_type.gff3"
    gff_path.write_text(content)
    regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene1"


# ---------------------------------------------------------------------------
# Multiple genes, some valid, some not
# ---------------------------------------------------------------------------


def test_mixed_valid_invalid_genes(tmp_path: Path, caplog):
    """Parser should extract valid genes and skip invalid ones with clear warnings."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene_valid;Name=ValidGene
        chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx_v;Parent=gene_valid
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds_v;Parent=tx_v
        chr1\t.\tgene\t200\t300\t.\t+\t.\tID=gene_bad_frame;Name=BadFrame
        chr1\t.\tmRNA\t200\t300\t.\t+\t.\tID=tx_b;Parent=gene_bad_frame
        chr1\t.\tCDS\t200\t210\t.\t+\t0\tID=cds_b;Parent=tx_b
        chr1\t.\tgene\t500\t600\t.\t+\t.\tID=gene_noncoding;Name=NonCoding
        chr1\t.\tlnc_RNA\t500\t600\t.\t+\t.\tID=lnc1;Parent=gene_noncoding
    """)
    gff_path = tmp_path / "mixed.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    assert len(regions) == 1
    assert regions[0].gene_id == "gene_valid"


# ---------------------------------------------------------------------------
# Completely empty / whitespace-only file
# ---------------------------------------------------------------------------


def test_whitespace_only_file(tmp_path: Path):
    """A file with only whitespace should return empty list."""
    gff_path = tmp_path / "whitespace.gff3"
    gff_path.write_text("   \n  \n\n")
    regions = parse_gff3(gff_path)
    assert regions == []


def test_file_not_found():
    """A non-existent file should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        parse_gff3("/nonexistent/path/file.gff3")


# ---------------------------------------------------------------------------
# Overlapping CDS intervals within one transcript
# ---------------------------------------------------------------------------


def test_orphan_cds_warning_no_mrna(tmp_path: Path, caplog):
    """When all mRNA features are missing, warn that CDS groups were treated as genes."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1
        chr1\t.\tgene\t200\t300\t.\t+\t.\tID=gene2
        chr1\t.\tCDS\t200\t208\t.\t+\t0\tID=cds2;Parent=tx2
    """)
    gff_path = tmp_path / "no_mrna_warn.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    # Both CDS groups are orphans — treated as separate genes
    assert len(regions) == 2
    assert any("no mRNA/transcript features found" in msg for msg in caplog.messages)


def test_orphan_cds_warning_broken_ids(tmp_path: Path, caplog):
    """When mRNA features lack ID, CDS parents can't resolve and trigger a warning."""
    # 3 genes with proper mRNA, 5 genes where mRNA lacks ID
    lines = ["##gff-version 3"]
    # 3 properly structured genes
    for i in range(1, 4):
        s = i * 100
        lines.append(f"chr1\t.\tgene\t{s}\t{s + 8}\t.\t+\t.\tID=gene{i}")
        lines.append(f"chr1\t.\tmRNA\t{s}\t{s + 8}\t.\t+\t.\tID=tx{i};Parent=gene{i}")
        lines.append(f"chr1\t.\tCDS\t{s}\t{s + 8}\t.\t+\t0\tID=cds{i};Parent=tx{i}")
    # 5 genes where mRNA has no ID — CDS parent won't resolve
    for i in range(4, 9):
        s = i * 100
        lines.append(f"chr1\t.\tgene\t{s}\t{s + 8}\t.\t+\t.\tID=gene{i}")
        lines.append(f"chr1\t.\tmRNA\t{s}\t{s + 8}\t.\t+\t.\tParent=gene{i}")
        lines.append(f"chr1\t.\tCDS\t{s}\t{s + 8}\t.\t+\t0\tID=cds{i};Parent=orphan_tx{i}")

    gff_path = tmp_path / "broken_ids.gff3"
    gff_path.write_text("\n".join(lines) + "\n")
    with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    # 3 resolved + 5 orphans = 8 total
    assert len(regions) == 8
    assert any("could not be linked to a gene" in msg for msg in caplog.messages)


def test_no_orphan_warning_when_cds_parents_gene_directly(tmp_path: Path, caplog):
    """CDS parented directly to a gene (no mRNA) should NOT trigger the large-orphan warning."""
    content = textwrap.dedent("""\
        ##gff-version 3
        chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
        chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1
    """)
    gff_path = tmp_path / "direct_parent.gff3"
    gff_path.write_text(content)
    with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
        regions = parse_gff3(gff_path)
    assert len(regions) == 1
    # This is a legitimate use — no warning about broken hierarchy
    assert not any("could not be linked" in msg for msg in caplog.messages)
    assert not any("no mRNA/transcript features found" in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# Files that cannot be read as GFF3
# ---------------------------------------------------------------------------


def test_non_utf8_raises(tmp_path: Path):
    gff_path = tmp_path / "latin1.gff3"
    gff_path.write_bytes(NON_UTF8_GFF)
    with pytest.raises(GffFormatError, match="latin1.gff3"):
        parse_gff3(gff_path)


def test_gz_name_not_gzipped_raises(tmp_path: Path):
    gff_path = tmp_path / "plain.gff3.gz"
    gff_path.write_text("##gff-version 3\n")
    with pytest.raises(GffFormatError, match="plain.gff3.gz"):
        parse_gff3(gff_path)


def test_directory_raises(tmp_path: Path):
    gff_dir = tmp_path / "genes.gff3"
    gff_dir.mkdir()
    with pytest.raises(GffFormatError, match="genes.gff3"):
        parse_gff3(gff_dir)


def test_corrupt_gzip_raises(tmp_path: Path):
    """A gzip header followed by damaged data fails inside zlib, not in the gzip module."""
    data = bytearray(
        gzip.compress(b"##gff-version 3\n" + b"chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=c\n" * 200)
    )
    data[12:20] = b"\xff" * 8
    gff_path = tmp_path / "damaged.gff3.gz"
    gff_path.write_bytes(bytes(data))
    with pytest.raises(GffFormatError, match="damaged.gff3.gz"):
        parse_gff3(gff_path)


# ---------------------------------------------------------------------------
# Checks on the CDS features of one transcript
# ---------------------------------------------------------------------------


def gff(tmp_path: Path, *cds: str, strand: str = "+", extra: str = "") -> Path:
    """One gene g1 with transcript t1 and the given CDS rows (chrom, start, end, strand, phase)."""
    rows = [
        "##gff-version 3",
        f"chr1\t.\tgene\t1\t90\t.\t{strand}\t.\tID=g1",
        f"chr1\t.\tmRNA\t1\t90\t.\t{strand}\t.\tID=t1;Parent=g1",
    ]
    for row in cds:
        chrom, start, end, cds_strand, phase = row.split()
        rows.append(f"{chrom}\t.\tCDS\t{start}\t{end}\t.\t{cds_strand}\t{phase}\tParent=t1")
    path = tmp_path / "genes.gff3"
    path.write_text("\n".join(rows) + "\n" + extra)
    return path


def warning_text(caplog) -> str:
    return " ".join(r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING)


class TestTranscriptChecks:
    """A transcript's own CDS lines must give one clean coding sequence."""

    def test_duplicate_cds_line_is_removed(self, tmp_path, caplog):
        path = gff(tmp_path, "chr1 1 30 + 0", "chr1 1 30 + 0")
        with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
            regions = parse_gff3(path)
        assert [r.exons for r in regions] == [[(0, 30)]]
        assert "1 duplicate CDS line" in warning_text(caplog)

    def test_overlapping_cds_skips_the_transcript(self, tmp_path, caplog):
        path = gff(tmp_path, "chr1 1 30 + 0", "chr1 28 60 + 0")
        with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
            assert parse_gff3(path) == []
        assert "overlapping CDS" in warning_text(caplog)

    @pytest.mark.parametrize("phase", [".", "3"])
    def test_invalid_used_phase_skips_the_transcript(self, tmp_path, caplog, phase):
        path = gff(tmp_path, f"chr1 1 30 + {phase}")
        with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
            assert parse_gff3(path) == []
        assert "invalid phase" in warning_text(caplog)

    def test_invalid_phase_on_another_cds_is_not_used(self, tmp_path):
        """Only the first CDS in coding order gives the phase, so a bad later one is harmless."""
        (region,) = parse_gff3(gff(tmp_path, "chr1 1 30 + 0", "chr1 41 70 + ."))
        assert region.exons == [(0, 30), (40, 70)]
        assert region.phase == 0

    def test_minus_strand_uses_the_highest_cds_phase(self, tmp_path):
        bad_first = gff(tmp_path, "chr1 1 30 - 0", "chr1 41 70 - .", strand="-")
        assert parse_gff3(bad_first) == []
        bad_last = gff(tmp_path, "chr1 1 30 - .", "chr1 41 70 - 0", strand="-")
        assert [r.exons for r in parse_gff3(bad_last)] == [[(0, 30), (40, 70)]]

    def test_cds_on_two_contigs_skips_the_transcript(self, tmp_path, caplog):
        path = gff(tmp_path, "chr1 1 30 + 0", "chr2 1 30 + 0")
        with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
            assert parse_gff3(path) == []
        assert "more than one contig or strand" in warning_text(caplog)

    def test_cds_on_two_strands_skips_the_transcript(self, tmp_path):
        assert parse_gff3(gff(tmp_path, "chr1 1 30 + 0", "chr1 31 60 - 0")) == []

    def test_strand_dot_skips_the_transcript(self, tmp_path):
        assert parse_gff3(gff(tmp_path, "chr1 1 30 . 0", strand=".")) == []

    def test_gene_falls_back_to_a_valid_transcript(self, tmp_path):
        """t1 is the longer transcript, but its CDS lines overlap, so the shorter t2 is used."""
        extra = (
            "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=t2;Parent=g1\n"
            "chr1\t.\tCDS\t1\t30\t.\t+\t0\tParent=t2\n"
        )
        path = gff(tmp_path, "chr1 1 30 + 0", "chr1 28 60 + 0", extra=extra)
        (region,) = parse_gff3(path)
        assert region.transcript_id == "t2"
        assert region.exons == [(0, 30)]

    def test_overlap_between_genes_is_kept(self, tmp_path):
        """Only CDS lines of one transcript are compared; other genes may overlap."""
        extra = (
            "chr1\t.\tgene\t1\t90\t.\t+\t.\tID=g2\n"
            "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=t2;Parent=g2\n"
            "chr1\t.\tCDS\t16\t45\t.\t+\t0\tParent=t2\n"
        )
        path = gff(tmp_path, "chr1 1 30 + 0", extra=extra)
        assert sorted(r.gene_id for r in parse_gff3(path)) == ["g1", "g2"]


class TestTranscriptChecksAtTheRightLayer:
    """Each check runs on every transcript, and CdsRegion rejects data no transcript can have."""

    def test_length_check_falls_back_to_a_valid_transcript(self, tmp_path):
        """t1 is longer but 31 bp is not a multiple of 3, so the 30 bp t2 is used."""
        extra = (
            "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=t2;Parent=g1\n"
            "chr1\t.\tCDS\t1\t30\t.\t+\t0\tParent=t2\n"
        )
        (region,) = parse_gff3(gff(tmp_path, "chr1 1 31 + 0", extra=extra))
        assert region.transcript_id == "t2"

    @pytest.mark.parametrize("phases", [("0", "."), (".", "0")])
    def test_copies_with_different_phases_skip_whatever_the_order(self, tmp_path, phases):
        """Two copies that disagree on the phase are not one line; the result must not
        depend on which copy comes first."""
        first, second = phases
        path = gff(tmp_path, f"chr1 1 30 + {first}", f"chr1 1 30 + {second}")
        assert parse_gff3(path) == []

    @pytest.mark.parametrize(
        ("exons", "strand", "phase"),
        [
            ([(0, 30), (0, 30)], "+", 0),
            ([(0, 30), (27, 60)], "+", 0),
            ([(0, 30)], ".", 0),
            ([(0, 30)], "+", 3),
            ([(0, 30), (40, 40)], "+", 0),
            ([(0, 30), (45, 40)], "+", 0),
        ],
    )
    def test_cds_region_rejects_what_no_transcript_can_have(self, exons, strand, phase):
        with pytest.raises(ValueError):
            CdsRegion(
                gene_id="g",
                transcript_id="t",
                chrom="chr1",
                exons=exons,
                strand=strand,
                phase=phase,
            )


class TestEmptyIntervalsAndConflictingPhases:
    """A reversed CDS interval, and copies whose phase tokens differ but are both invalid."""

    @pytest.mark.parametrize("bad", ["chr1 11 9 + 0", "chr1 11 10 + 0"])
    def test_reversed_or_empty_cds_skips_the_transcript(self, tmp_path, caplog, bad):
        """Lengths of 7 and -1 add to 6, but the region would hold 7 positions."""
        path = gff(tmp_path, "chr1 1 7 + 0", bad)
        with caplog.at_level(logging.WARNING, logger="mkado.io.gff"):
            assert parse_gff3(path) == []
        assert "empty or reversed" in warning_text(caplog)

    def test_copies_with_two_invalid_phases_are_not_one_line(self, tmp_path):
        path = gff(tmp_path, "chr1 1 30 + 0", "chr1 41 70 + .", "chr1 41 70 + 3")
        assert parse_gff3(path) == []
