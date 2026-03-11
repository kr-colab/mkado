"""Tests for CDS coordinate mapping."""

from __future__ import annotations

import pytest

from mkado.core.cds import CdsRegion


class TestCdsRegionPlusStrand:
    """Tests for plus-strand CDS regions."""

    def test_simple_single_exon(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 109)],  # 9 bases = 3 codons
            strand="+",
        )
        assert cds.cds_length == 9
        assert cds.num_codons() == 3
        assert cds.is_valid()

    def test_genomic_positions_plus_strand(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 106)],  # 6 bases = 2 codons
            strand="+",
        )
        positions = cds.genomic_positions()
        assert positions == [100, 101, 102, 103, 104, 105]

    def test_codon_positions(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 109)],
            strand="+",
        )
        assert cds.codon_positions(0) == (100, 101, 102)
        assert cds.codon_positions(1) == (103, 104, 105)
        assert cds.codon_positions(2) == (106, 107, 108)

    def test_genomic_pos_to_codon_index(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 109)],
            strand="+",
        )
        assert cds.genomic_pos_to_codon_index(100) == 0
        assert cds.genomic_pos_to_codon_index(102) == 0
        assert cds.genomic_pos_to_codon_index(103) == 1
        assert cds.genomic_pos_to_codon_index(108) == 2
        assert cds.genomic_pos_to_codon_index(99) is None
        assert cds.genomic_pos_to_codon_index(109) is None

    def test_genomic_pos_to_codon_offset(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 109)],
            strand="+",
        )
        assert cds.genomic_pos_to_codon_offset(100) == 0
        assert cds.genomic_pos_to_codon_offset(101) == 1
        assert cds.genomic_pos_to_codon_offset(102) == 2
        assert cds.genomic_pos_to_codon_offset(103) == 0

    def test_contains_position(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 106)],
            strand="+",
        )
        assert cds.contains_position(100)
        assert cds.contains_position(105)
        assert not cds.contains_position(99)
        assert not cds.contains_position(106)

    def test_multi_exon_plus_strand(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 104), (200, 205)],  # 4 + 5 = 9 bases = 3 codons
            strand="+",
        )
        assert cds.cds_length == 9
        assert cds.num_codons() == 3
        positions = cds.genomic_positions()
        assert positions == [100, 101, 102, 103, 200, 201, 202, 203, 204]
        # First codon spans exon boundary
        assert cds.codon_positions(0) == (100, 101, 102)
        # Second codon: last base of exon1 + first two of exon2
        assert cds.codon_positions(1) == (103, 200, 201)

    def test_invalid_cds_length(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 108)],  # 8 bases, not divisible by 3
            strand="+",
        )
        assert not cds.is_valid()


class TestCdsRegionMinusStrand:
    """Tests for minus-strand CDS regions."""

    def test_minus_strand_positions_reversed(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 106)],  # 6 bases = 2 codons
            strand="-",
        )
        positions = cds.genomic_positions()
        # Minus strand: positions are reversed
        assert positions == [105, 104, 103, 102, 101, 100]

    def test_minus_strand_codon_positions(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 106)],
            strand="-",
        )
        # First codon on minus strand: highest genomic positions
        assert cds.codon_positions(0) == (105, 104, 103)
        assert cds.codon_positions(1) == (102, 101, 100)

    def test_minus_strand_codon_index_mapping(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 106)],
            strand="-",
        )
        # Position 105 is first base of first codon (minus strand)
        assert cds.genomic_pos_to_codon_index(105) == 0
        assert cds.genomic_pos_to_codon_index(103) == 0
        assert cds.genomic_pos_to_codon_index(102) == 1
        assert cds.genomic_pos_to_codon_index(100) == 1

    def test_minus_strand_multi_exon(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 103), (200, 206)],  # 3 + 6 = 9 bases
            strand="-",
        )
        assert cds.num_codons() == 3
        positions = cds.genomic_positions()
        # Reversed: exon2 positions first (higher), then exon1
        assert positions == [205, 204, 203, 202, 201, 200, 102, 101, 100]


class TestCdsRegionPhase:
    """Tests for CDS regions with non-zero phase."""

    def test_phase_1_skips_first_base(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 107)],  # 7 bases, phase 1 -> skip 1 -> 6 bases = 2 codons
            strand="+",
            phase=1,
        )
        assert cds.cds_length == 6
        assert cds.num_codons() == 2
        positions = cds.genomic_positions()
        assert positions == [101, 102, 103, 104, 105, 106]

    def test_phase_2_skips_two_bases(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(100, 108)],  # 8 bases, phase 2 -> skip 2 -> 6 bases = 2 codons
            strand="+",
            phase=2,
        )
        assert cds.cds_length == 6
        assert cds.num_codons() == 2
        positions = cds.genomic_positions()
        assert positions == [102, 103, 104, 105, 106, 107]


class TestCdsRegionExtractCodon:
    """Tests for codon extraction with a mock ref_fetch."""

    def test_extract_codon_plus_strand(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(0, 9)],
            strand="+",
        )
        # Simulate reference: ATGCCCAAA
        ref_seq = "ATGCCCAAA"

        def ref_fetch(chrom, pos):
            return ref_seq[pos]

        assert cds.extract_codon(0, ref_fetch) == "ATG"
        assert cds.extract_codon(1, ref_fetch) == "CCC"
        assert cds.extract_codon(2, ref_fetch) == "AAA"

    def test_extract_codon_minus_strand(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(0, 6)],
            strand="-",
        )
        # Reference forward: AACCGG (indices 0-5)
        # Minus strand reads from pos 5 backwards: G(5) G(4) C(3) C(2) A(1) A(0)
        # Complement:                              C    C    G    G    T    T
        ref_seq = "AACCGG"

        def ref_fetch(chrom, pos):
            return ref_seq[pos]

        # First codon (positions 5,4,3): G,G,C -> complement: C,C,G
        assert cds.extract_codon(0, ref_fetch) == "CCG"
        # Second codon (positions 2,1,0): C,A,A -> complement: G,T,T
        assert cds.extract_codon(1, ref_fetch) == "GTT"


class TestExonOrdering:
    """Tests that exon ordering is handled correctly."""

    def test_unsorted_exons_are_sorted(self):
        cds = CdsRegion(
            gene_id="gene1",
            transcript_id="tx1",
            chrom="chr1",
            exons=[(200, 203), (100, 103)],  # Given in wrong order
            strand="+",
        )
        # Should sort to [(100, 103), (200, 203)]
        assert cds.exons == [(100, 103), (200, 203)]
        positions = cds.genomic_positions()
        assert positions == [100, 101, 102, 200, 201, 202]
