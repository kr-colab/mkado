"""CDS coordinate mapping for VCF-based analysis."""

from __future__ import annotations

from dataclasses import dataclass, field

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass
class CdsRegion:
    """A coding sequence region defined by exon coordinates.

    Coordinates are 0-based half-open (like BED/Python).
    For minus-strand genes, exons should be provided in genomic order;
    they will be processed in reverse order internally.
    """

    gene_id: str
    transcript_id: str
    chrom: str
    exons: list[tuple[int, int]]
    """CDS exon intervals, 0-based half-open, in genomic order."""
    strand: str
    """'+' or '-'."""
    phase: int = 0
    """Phase of the first CDS exon (0, 1, or 2)."""

    _positions: list[int] = field(default=None, init=False, repr=False)
    _pos_to_idx: dict[int, int] = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        self.exons = sorted(self.exons, key=lambda e: e[0])

    def _build_positions(self) -> None:
        if self._positions is not None:
            return
        positions: list[int] = []
        for start, end in self.exons:
            positions.extend(range(start, end))
        if self.strand == "-":
            positions = positions[::-1]
        # Skip phase offset bases from the beginning
        if self.phase > 0:
            positions = positions[self.phase :]
        self._positions = positions
        self._pos_to_idx = {pos: i for i, pos in enumerate(positions)}

    @property
    def cds_length(self) -> int:
        """Total CDS length in bases (after phase adjustment)."""
        self._build_positions()
        return len(self._positions)

    def num_codons(self) -> int:
        """Number of complete codons."""
        return self.cds_length // 3

    def is_valid(self) -> bool:
        """Check if CDS length is divisible by 3."""
        return self.cds_length % 3 == 0

    def genomic_positions(self) -> list[int]:
        """All CDS positions in codon reading order.

        For plus-strand genes, this is ascending genomic order.
        For minus-strand genes, this is descending genomic order.
        """
        self._build_positions()
        return list(self._positions)

    def codon_positions(self, codon_index: int) -> tuple[int, int, int]:
        """Get the 3 genomic positions for a codon.

        Args:
            codon_index: 0-based codon index in the CDS.

        Returns:
            Tuple of 3 genomic positions (0-based) for the codon.
        """
        self._build_positions()
        base = codon_index * 3
        return (self._positions[base], self._positions[base + 1], self._positions[base + 2])

    def genomic_pos_to_codon_index(self, pos: int) -> int | None:
        """Map a genomic position to its codon index.

        Args:
            pos: 0-based genomic position.

        Returns:
            Codon index, or None if position is not in this CDS.
        """
        self._build_positions()
        idx = self._pos_to_idx.get(pos)
        if idx is None:
            return None
        return idx // 3

    def genomic_pos_to_codon_offset(self, pos: int) -> int | None:
        """Map a genomic position to its offset within the codon (0, 1, or 2).

        Args:
            pos: 0-based genomic position.

        Returns:
            Offset within the codon, or None if position is not in this CDS.
        """
        self._build_positions()
        idx = self._pos_to_idx.get(pos)
        if idx is None:
            return None
        return idx % 3

    def contains_position(self, pos: int) -> bool:
        """Check if a genomic position falls within this CDS."""
        self._build_positions()
        return pos in self._pos_to_idx

    def extract_codon(self, codon_index: int, ref_fetch: callable) -> str:
        """Extract a reference codon sequence.

        Args:
            codon_index: 0-based codon index.
            ref_fetch: Callable that takes (chrom, pos) and returns a base.

        Returns:
            3-letter codon string (uppercase). For minus strand, reverse-complemented.
        """
        p1, p2, p3 = self.codon_positions(codon_index)
        bases = ref_fetch(self.chrom, p1) + ref_fetch(self.chrom, p2) + ref_fetch(self.chrom, p3)
        if self.strand == "-":
            bases = bases.translate(_COMPLEMENT)
        return bases.upper()
