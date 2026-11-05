"""FASTA file parser and writer."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

from mkado.errors import InputError


class FastaFormatError(InputError):
    """A FASTA file has a line that cannot be read as FASTA."""


# The digit check runs on every sequence line, so it is a compiled search that runs in C, not
# a Python loop over each base.
_DIGIT = re.compile(r"\d")


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Parse a FASTA file and yield (name, sequence) tuples.

    Whitespace inside a sequence line is removed, so a sequence written in blocks of
    bases reads the same as one with no spaces.

    Args:
        path: Path to FASTA file

    Yields:
        Tuples of (sequence_name, sequence_string)

    Raises:
        FastaFormatError: A sequence line has a digit or comes before the first header,
            or a header has no name.
    """
    path = Path(path)
    name: str | None = None
    seq_parts: list[str] = []

    def error(line_number: int, problem: str) -> FastaFormatError:
        return FastaFormatError(f"{path.name}, line {line_number}: {problem}")

    with path.open() as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts)
                words = line[1:].split()
                if not words:
                    raise error(line_number, "the header has no name.")
                name = words[0]
                seq_parts = []
                continue
            if name is None:
                raise error(line_number, "sequence data comes before the first '>' header.")
            # Each whitespace character would otherwise count as a base and shift every
            # codon after it.
            bases = "".join(line.split())
            # A digit is not missing data: it marks a numbered layout such as GenBank, and
            # reading it as a base would shift the codons in the same way.
            if _DIGIT.search(bases):
                raise error(
                    line_number,
                    f"record {name} has a digit in its sequence. Numbered sequence lines, "
                    "as in GenBank files, are not FASTA.",
                )
            seq_parts.append(bases.upper())

        if name is not None:
            yield name, "".join(seq_parts)


def write_fasta(
    sequences: list[tuple[str, str]],
    path: str | Path,
    line_width: int = 70,
) -> None:
    """Write sequences to a FASTA file.

    Args:
        sequences: List of (name, sequence) tuples
        path: Output file path
        line_width: Characters per line for sequence wrapping
    """
    path = Path(path)
    with path.open("w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")
