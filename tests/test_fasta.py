"""Tests for FASTA parser."""

from pathlib import Path

import pytest

from mkado.core.sequences import AlignmentLengthError
from mkado.errors import InputError
from mkado.io.fasta import FastaFormatError, read_fasta, write_fasta
from mkado.io.gff import GffFormatError


def test_read_fasta(tmp_path: Path) -> None:
    """Test reading a FASTA file."""
    fasta_content = """>seq1
ATGCATGC
ATGCATGC
>seq2
GGGGCCCC
"""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    sequences = list(read_fasta(fasta_file))

    assert len(sequences) == 2
    assert sequences[0] == ("seq1", "ATGCATGCATGCATGC")
    assert sequences[1] == ("seq2", "GGGGCCCC")


def test_read_fasta_with_header_info(tmp_path: Path) -> None:
    """Test that only the first word of the header is used as name."""
    fasta_content = """>seq1 description here
ATGCATGC
"""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    sequences = list(read_fasta(fasta_file))

    assert sequences[0][0] == "seq1"


def test_read_fasta_uppercase(tmp_path: Path) -> None:
    """Test that sequences are converted to uppercase."""
    fasta_content = """>seq1
atgcATGC
"""
    fasta_file = tmp_path / "test.fa"
    fasta_file.write_text(fasta_content)

    sequences = list(read_fasta(fasta_file))

    assert sequences[0][1] == "ATGCATGC"


def test_write_fasta(tmp_path: Path) -> None:
    """Test writing a FASTA file."""
    sequences = [("seq1", "ATGCATGC"), ("seq2", "GGGGCCCC")]
    fasta_file = tmp_path / "output.fa"

    write_fasta(sequences, fasta_file, line_width=4)

    content = fasta_file.read_text()
    lines = content.strip().split("\n")

    assert lines[0] == ">seq1"
    assert lines[1] == "ATGC"
    assert lines[2] == "ATGC"
    assert lines[3] == ">seq2"


def test_read_empty_fasta(tmp_path: Path) -> None:
    """Test reading an empty FASTA file."""
    fasta_file = tmp_path / "empty.fa"
    fasta_file.write_text("")

    sequences = list(read_fasta(fasta_file))

    assert len(sequences) == 0


def read(tmp_path: Path, text: str, name: str = "g.fa") -> list[tuple[str, str]]:
    path = tmp_path / name
    path.write_bytes(text.encode())
    return list(read_fasta(path))


@pytest.mark.parametrize(
    ("line", "bases"),
    [
        ("ATG AAA", "ATGAAA"),
        ("ATG\tAAA", "ATGAAA"),
        # Some tools write a sequence in blocks of 10 bases with a space between the blocks.
        ("ATGAAACCCG GGTTTAAACC", "ATGAAACCCGGGTTTAAACC"),
    ],
)
def test_whitespace_inside_sequence_is_removed(tmp_path: Path, line: str, bases: str) -> None:
    assert read(tmp_path, f">a\n{line}\n") == [("a", bases)]


def test_windows_line_endings(tmp_path: Path) -> None:
    assert read(tmp_path, ">a\r\nATGAAA\r\nCCC\r\n") == [("a", "ATGAAACCC")]


def test_digit_raises(tmp_path: Path) -> None:
    """A numbered sequence line is not FASTA; reading the digits as bases would shift codons."""
    with pytest.raises(FastaFormatError) as excinfo:
        read(tmp_path, ">a\n1 atgaaa\n", name="numbered.fa")
    message = str(excinfo.value)
    assert "numbered.fa" in message
    assert "line 2" in message
    assert "record a" in message


def test_sequence_before_first_header_raises(tmp_path: Path) -> None:
    with pytest.raises(FastaFormatError, match="line 1"):
        read(tmp_path, "ATGAAA\n>a\nATGAAA\n")


@pytest.mark.parametrize("header", [">", ">   "])
def test_header_without_name_raises(tmp_path: Path, header: str) -> None:
    with pytest.raises(FastaFormatError, match="line 3"):
        read(tmp_path, f">a\nATGAAA\n{header}\nATGAAG\n")


def test_input_errors_are_value_errors() -> None:
    """Both input errors share one base class, and code that catches ValueError still works."""
    for error in (FastaFormatError, AlignmentLengthError, GffFormatError):
        assert issubclass(error, InputError)
        assert issubclass(error, ValueError)


def test_non_utf8_raises(tmp_path: Path) -> None:
    path = tmp_path / "latin1.fa"
    path.write_bytes(b">a\xff\nATGAAA\n")
    with pytest.raises(FastaFormatError, match="latin1.fa"):
        list(read_fasta(path))


def test_directory_raises(tmp_path: Path) -> None:
    fasta_dir = tmp_path / "genes.fa"
    fasta_dir.mkdir()
    with pytest.raises(FastaFormatError, match="genes.fa"):
        list(read_fasta(fasta_dir))
