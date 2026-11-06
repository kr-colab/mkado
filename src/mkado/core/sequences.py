"""Sequence and SequenceSet classes for handling aligned sequences."""

from __future__ import annotations

import numbers
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.data.genetic_codes import CLEAN_BASES, CLEAN_CODONS
from mkado.errors import InputError
from mkado.io.fasta import read_fasta

if TYPE_CHECKING:
    pass


# A base that counts as data, used to find whether a sequence goes on after a stop codon.
_BASE = re.compile("[ACGT]")


class AlignmentLengthError(InputError):
    """Sequences that must be aligned have different lengths."""


# A file can hold thousands of sequences, and an unaligned file can have as many lengths,
# so the error lists only this many lengths, and this many names for each length.
_LENGTHS_SHOWN = 5
_NAMES_PER_LENGTH = 3


def _first_few(shown: list[str], total: int, noun: str) -> str:
    """``shown`` joined with commas, then a count of the ``total`` items not shown."""
    text = ", ".join(shown)
    if total > len(shown):
        text += f" and {total - len(shown)} more {noun}"
    return text


def _length_mismatch_message(sequences: list[Sequence]) -> str:
    """The lengths in ``sequences``, most common first, with the first few names that have each."""
    names_by_length: dict[int, list[str]] = {}
    for seq in sequences:
        names_by_length.setdefault(len(seq), []).append(seq.name)
    # The most common length is usually the alignment length, so it comes first and the
    # sequences that differ from it follow. A tie goes to the shorter length.
    ordered = sorted(names_by_length.items(), key=lambda item: (-len(item[1]), item[0]))
    groups = [
        f"{length} bp ({_first_few(names[:_NAMES_PER_LENGTH], len(names), 'names')})"
        for length, names in ordered[:_LENGTHS_SHOWN]
    ]
    return "Sequences have different lengths: " + _first_few(groups, len(ordered), "lengths")


def as_sequence_set(
    source: SequenceSet | str | Path,
    reading_frame: int | None,
    genetic_code: GeneticCode,
) -> SequenceSet:
    """The SequenceSet for an analysis argument that is a set or a FASTA path.

    A path is read in ``reading_frame``, or in frame 1 when it is None. A set already has
    a frame, so a ``reading_frame`` that differs from it raises ValueError: the set would
    otherwise be read in its own frame and the argument would have no effect.
    """
    if isinstance(source, (str, Path)):
        frame = 1 if reading_frame is None else reading_frame
        return SequenceSet.from_fasta(source, frame, genetic_code)
    if reading_frame is not None and reading_frame != source.reading_frame:
        raise ValueError(
            f"reading_frame={reading_frame} was given, but the SequenceSet is in frame "
            f"{source.reading_frame}"
        )
    return source


def position_bases(counts: dict[str, int], position: int) -> dict[str, int]:
    """Copies of each base at one codon position, from codon counts."""
    bases: dict[str, int] = {}
    for codon, count in counts.items():
        bases[codon[position]] = bases.get(codon[position], 0) + count
    return bases


def ancestral_base(bases: dict[str, int], source_codons: set[str], position: int) -> str | None:
    """The base shared with ``source_codons`` at ``position`` that is most common in ``bases``.

    A tie goes to the first base in alphabetical order, so the result does not depend
    on the order of the sequences. None when no base is shared.
    """
    shared = bases.keys() & {codon[position] for codon in source_codons}
    if not shared:
        return None
    return min(shared, key=lambda base: (-bases[base], base))


@dataclass
class Sequence:
    """Represents a single biological sequence."""

    name: str
    sequence: str

    def __post_init__(self) -> None:
        # The clean-codon check matches uppercase bases only, so lowercase input
        # (soft-masked sequence, for example) would otherwise all be missing data.
        self.sequence = self.sequence.upper()

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: int | slice) -> str:
        return self.sequence[index]

    def get_codon(self, codon_index: int, reading_frame: int = 1) -> str:
        """Get the codon at a given index.

        Args:
            codon_index: Zero-based codon index
            reading_frame: Reading frame (1, 2, or 3)

        Returns:
            Three-letter codon string
        """
        start = (reading_frame - 1) + (codon_index * 3)
        return self.sequence[start : start + 3]

    def num_codons(self, reading_frame: int = 1) -> int:
        """Get the number of complete codons in the sequence."""
        return (len(self.sequence) - (reading_frame - 1)) // 3


@dataclass
class SequenceSet:
    """Represents a set of aligned sequences (e.g., from one species)."""

    sequences: list[Sequence] = field(default_factory=list)
    reading_frame: int = 1
    genetic_code: GeneticCode = field(default_factory=lambda: DEFAULT_CODE)
    _codon_set_clean_cache: dict[int, set[str]] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )
    """Per-position memoization for ``codon_set_clean()``. Treat the SequenceSet
    as immutable once any caching method has been called; mutating ``sequences``
    or ``reading_frame`` afterwards will leave stale entries here."""

    def __post_init__(self) -> None:
        # Any other value reads codons at the wrong columns or fails later: -1 slices from
        # the end, and a bool or a float passes a plain comparison with 1, 2 or 3.
        frame = self.reading_frame
        if (
            isinstance(frame, bool)
            or not isinstance(frame, numbers.Integral)
            or frame not in (1, 2, 3)
        ):
            raise ValueError(f"reading_frame must be 1, 2 or 3, not {frame!r}")
        # Every length is read from the first sequence, so a set of mixed lengths would
        # give results that depend on which sequence comes first.
        if len({len(seq) for seq in self.sequences}) > 1:
            raise AlignmentLengthError(_length_mismatch_message(self.sequences))

    @classmethod
    def from_fasta(
        cls,
        path: str | Path,
        reading_frame: int = 1,
        genetic_code: GeneticCode | None = None,
    ) -> SequenceSet:
        """Load sequences from a FASTA file.

        Args:
            path: Path to FASTA file
            reading_frame: Reading frame (1, 2, or 3)
            genetic_code: Genetic code for translation

        Returns:
            SequenceSet containing all sequences from the file
        """
        sequences = [Sequence(name, seq) for name, seq in read_fasta(path)]
        try:
            return cls(
                sequences=sequences,
                reading_frame=reading_frame,
                genetic_code=genetic_code or DEFAULT_CODE,
            )
        except AlignmentLengthError as e:
            raise AlignmentLengthError(f"{Path(path).name}: {e}") from None

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> Sequence:
        return self.sequences[index]

    def __getstate__(self) -> dict:
        # Drop the per-position cache when pickling — it can be re-populated
        # cheaply in the receiver and would otherwise inflate pickle size by
        # O(num_codons) set objects.
        state = self.__dict__.copy()
        state["_codon_set_clean_cache"] = {}
        return state

    @property
    def num_codons(self) -> int:
        """Number of complete codons in each sequence."""
        if not self.sequences:
            return 0
        return self.sequences[0].num_codons(self.reading_frame)

    @property
    def alignment_length(self) -> int:
        """Length of the alignment in bases."""
        if not self.sequences:
            return 0
        return len(self.sequences[0])

    def codon_set(self, codon_index: int) -> set[str]:
        """Get the set of unique codons at a position across all sequences.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique codon strings
        """
        codons = set()
        for seq in self.sequences:
            codon = seq.get_codon(codon_index, self.reading_frame)
            codons.add(codon)
        return codons

    def codon_set_clean(self, codon_index: int) -> set[str]:
        """Get unique codons at a position, excluding any with a base other than A, C, G or T.

        Result is cached per ``codon_index`` for the lifetime of the
        ``SequenceSet``. Callers must not mutate the returned set.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique valid codon strings
        """
        cached = self._codon_set_clean_cache.get(codon_index)
        if cached is not None:
            return cached
        codons = self.codon_set(codon_index)
        clean = codons & CLEAN_CODONS
        self._codon_set_clean_cache[codon_index] = clean
        return clean

    def amino_set(self, codon_index: int) -> set[str]:
        """Get the set of unique amino acids at a codon position.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique amino acid characters
        """
        codons = self.codon_set(codon_index)
        return {self.genetic_code.translate(c) for c in codons}

    def amino_set_clean(self, codon_index: int) -> set[str]:
        """Get unique amino acids at a position, excluding ambiguous codons.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Set of unique amino acid characters
        """
        codons = self.codon_set_clean(codon_index)
        return {self.genetic_code.translate(c) for c in codons}

    def is_polymorphic(self, codon_index: int) -> bool:
        """Check if a codon position is polymorphic (has >1 unique codon).

        Args:
            codon_index: Zero-based codon index

        Returns:
            True if polymorphic
        """
        codons = self.codon_set_clean(codon_index)
        return len(codons) > 1

    def site_set(self, site_index: int) -> set[str]:
        """Get the set of unique nucleotides at a site across all sequences.

        Args:
            site_index: Zero-based nucleotide position

        Returns:
            Set of unique nucleotide characters
        """
        return {seq.sequence[site_index] for seq in self.sequences if site_index < len(seq)}

    def site_set_clean(self, site_index: int) -> set[str]:
        """Get unique nucleotides at a site, excluding any base other than A, C, G or T.

        Args:
            site_index: Zero-based nucleotide position

        Returns:
            Set of unique valid nucleotide characters
        """
        sites = self.site_set(site_index)
        return sites & CLEAN_BASES

    def codon_array(self) -> np.ndarray:
        """Get a 2D array of codons (n_sequences x n_codons).

        Returns:
            NumPy array of codon strings
        """
        n_seqs = len(self.sequences)
        n_codons = self.num_codons
        arr = np.empty((n_seqs, n_codons), dtype="U3")
        for i, seq in enumerate(self.sequences):
            for j in range(n_codons):
                arr[i, j] = seq.get_codon(j, self.reading_frame)
        return arr

    def polymorphic_codons(self) -> list[int]:
        """Get indices of polymorphic codon positions.

        Returns:
            List of codon indices that are polymorphic
        """
        return [i for i in range(self.num_codons) if self.is_polymorphic(i)]

    def site_codon_counts(self, codon_index: int) -> tuple[dict[str, int], int]:
        """Count each codon at a position, over sequences whose codon is all A, C, G or T.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Tuple of (codon counts, number of sequences counted). The count is the
            denominator of the frequencies, which is below the sequence total wherever
            a sequence is missing data at this codon.
        """
        counts: dict[str, int] = {}
        total = 0
        for seq in self.sequences:
            codon = seq.get_codon(codon_index, self.reading_frame)
            if codon in CLEAN_CODONS:
                counts[codon] = counts.get(codon, 0) + 1
                total += 1
        return counts, total

    def site_frequency_spectrum(self, codon_index: int) -> dict[str, float]:
        """Calculate allele frequencies at a codon position.

        Args:
            codon_index: Zero-based codon index

        Returns:
            Dict mapping codon strings to their frequencies (0-1)
        """
        counts, total = self.site_codon_counts(codon_index)
        if total == 0:
            return {}
        return {codon: count / total for codon, count in counts.items()}

    def derived_state(
        self, ancestral_source: SequenceSet, codon_index: int
    ) -> tuple[float, int] | None:
        """Derive the frequency and copy count of the non-ancestral alleles at a codon.

        The ancestral codon is the one this set shares with ``ancestral_source`` at the
        highest frequency here. Both figures cover every non-ancestral codon together,
        so a codon carrying two different derived alleles reports a count of two.

        Args:
            ancestral_source: Sequences used to infer the ancestral codon. This is the
                outgroup for the standard test, and the second outgroup when polarizing.
            codon_index: Zero-based codon index.

        Returns:
            Tuple of (derived frequency, derived count), or None when the ancestral
            state cannot be determined.
        """
        counts, total = self.site_codon_counts(codon_index)
        if total == 0:
            return None

        shared = set(counts) & ancestral_source.codon_set_clean(codon_index)
        if not shared:
            return None

        ancestral = max(shared, key=lambda codon: counts[codon])
        derived_count = total - counts[ancestral]
        # Divide the count rather than subtracting the ancestral frequency, so a site
        # sits exactly on its own frequency and a threshold comparison is not off by a
        # rounding.
        return derived_count / total, derived_count

    def position_derived_state(
        self, ancestral_source: SequenceSet, codon_index: int, position: int
    ) -> tuple[float, int] | None:
        """Derive the frequency and copy count of the non-ancestral bases at one codon position.

        Polarizing one position, not the whole codon, means that a difference at another
        position of the codon does not stop the ancestral state from being found. The
        ancestral base is the one this set shares with ``ancestral_source`` at that
        position that is most common here; a tie goes to the first base in alphabetical
        order, so the result does not depend on the order of the sequences.

        Args:
            ancestral_source: Sequences used to infer the ancestral base.
            codon_index: Zero-based codon index.
            position: Position in the codon, 0, 1 or 2.

        Returns:
            Tuple of (derived frequency, derived count), or None when no base at that
            position is shared, including when ``ancestral_source`` has no clean codon.
        """
        counts, total = self.site_codon_counts(codon_index)
        bases = position_bases(counts, position)
        ancestral = ancestral_base(bases, ancestral_source.codon_set_clean(codon_index), position)
        if ancestral is None:
            return None
        derived_count = total - bases[ancestral]
        return derived_count / total, derived_count

    def with_stops_masked(self, genetic_code: GeneticCode) -> tuple[SequenceSet, set[int]]:
        """A copy with each stop codon of ``genetic_code`` as NNN, and the internal stop codons.

        NNN is missing data wherever a codon is counted, so a stop codon then adds nothing
        to divergence, polymorphism or the site counts. The length does not change.

        A stop codon is internal when its own sequence has a later codon with an A, C, G or
        T base. A stop that only gaps or N follow is that sequence's terminal stop, even when
        another sequence goes on.
        """
        n_codons = self.num_codons
        if not genetic_code.stop_codons or not self.sequences or n_codons == 0:
            return self, set()
        start = self.reading_frame - 1
        # One array pass over all sequences, because the set is masked for every pair that
        # an analysis builds. A character that is not ASCII becomes one "?" byte, so the
        # positions do not move.
        raw = "".join(seq.sequence for seq in self.sequences).encode("ascii", "replace")
        frame = np.frombuffer(raw, dtype=np.uint8).reshape(len(self.sequences), -1)
        codon_bytes = frame[:, start : start + 3 * n_codons].reshape(len(self.sequences), -1, 3)
        codons = np.ascontiguousarray(codon_bytes).view("S3")[..., 0]
        is_stop = np.isin(codons, [codon.encode() for codon in genetic_code.stop_codons])
        if not is_stop.any():
            return self, set()

        sequences = list(self.sequences)
        internal: set[int] = set()
        end = start + 3 * n_codons
        rows, stop_codons = np.nonzero(is_stop)
        # Most sequences carry one stop, the terminal one, so slicing around each stop is
        # cheaper than rebuilding the sequence from a list of bases. np.nonzero gives the
        # stops of a row in order, so a later stop is not masked yet when an earlier one is
        # tested, and its bases make the earlier stop internal.
        for row, codon in zip(rows.tolist(), stop_codons.tolist(), strict=True):
            text = sequences[row].sequence
            pos = start + 3 * codon
            if _BASE.search(text, pos + 3, end):
                internal.add(codon)
            sequences[row] = Sequence(sequences[row].name, text[:pos] + "NNN" + text[pos + 3 :])

        masked = SequenceSet(
            sequences=sequences, reading_frame=self.reading_frame, genetic_code=self.genetic_code
        )
        return masked, internal

    def filter_by_name(self, pattern: str) -> SequenceSet:
        """Filter sequences by name pattern (substring match).

        Args:
            pattern: Substring to match in sequence names

        Returns:
            New SequenceSet containing only matching sequences
        """
        filtered = [seq for seq in self.sequences if pattern in seq.name]
        # The new SequenceSet starts with an empty cache: the parent's cached
        # codon sets contain codons from sequences that aren't in the filtered
        # subset, so reusing them would yield wrong results.
        return SequenceSet(
            sequences=filtered,
            reading_frame=self.reading_frame,
            genetic_code=self.genetic_code,
        )
