"""Builders for the hand-written alignments used across the test suite."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import Sequence, SequenceSet


def sequence_set(
    seqs: list[str], reading_frame: int = 1, genetic_code: GeneticCode = DEFAULT_CODE
) -> SequenceSet:
    """Build a SequenceSet from sequence strings, named s0, s1, and so on."""
    return SequenceSet(
        sequences=[Sequence(name=f"s{i}", sequence=s) for i, s in enumerate(seqs)],
        reading_frame=reading_frame,
        genetic_code=genetic_code,
    )


# One CDS line with GTF attributes (key "value"), which the GFF3 reader rejects.
GTF_LINE = 'chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1";\n'
# A GFF3 file with one byte that is not UTF-8.
NON_UTF8_GFF = b"##gff-version 3\nchr1\t.\tCDS\t1\t9\t.\t+\t0\tID=c;Name=\xff\n"


def assert_clean_error(result, *fragments: str) -> None:
    """The command stopped with one error line that has every fragment, not a traceback."""
    assert result.exit_code == 1, result.output
    # An uncaught ValueError would mean the user saw a traceback.
    assert isinstance(result.exception, SystemExit), result.exception
    assert "Error:" in result.output
    for fragment in fragments:
        assert fragment in result.output


def write_alignment_file(path: Path, records: dict[str, str]) -> Path:
    """Write a FASTA alignment from a name-to-sequence mapping."""
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
    return path


# Two genes of six codons for the batch command. Four ingroup sequences each
# carry at most one singleton polymorphism, one outgroup carries the fixed
# differences, and a second outgroup matches the ingroup, so polarization puts
# every fixed difference on the outgroup lineage.
# Standard counts: g1 Dn=1 Ds=1 Pn=0 Ps=2, g2 Dn=0 Ds=1 Pn=1 Ps=1.
TWO_GENES: dict[str, dict[str, list[str]]] = {
    "g1": {
        "ingroup": [
            "ATGAAACCCGGGTTTAAA",
            "ATGAAGCCCGGGTTTAAA",
            "ATGAAACCCGGGTTCAAA",
            "ATGAAACCCGGGTTTAAA",
        ],
        "outgroup": ["ATGAAACCAGGGTTTAAC"],
        "outgroup2": ["ATGAAACCCGGGTTTAAA"],
    },
    "g2": {
        "ingroup": [
            "ATGCCCAAAGGGTTTAAA",
            "ATGCCCAAGGGGTTTAAA",
            "ATGCCCAAAGGGTTTAAA",
            "ATGCCCAAAGGGTTAAAA",
        ],
        "outgroup": ["ATGCCCAAAGGCTTTAAA"],
        "outgroup2": ["ATGCCCAAAGGGTTTAAA"],
    },
}

# Sequence-name tags for the combined layout, one per group.
SPECIES = {"ingroup": "speciesA", "outgroup": "speciesB", "outgroup2": "speciesC"}


@dataclass
class TwoGeneDirs:
    """The two-gene set written in both batch layouts."""

    combined: Path
    """One file per gene, sequences named <gene>_<species>_<n>."""
    separate: Path
    """Per gene: <gene>_ingroup.fa, <gene>_outgroup.fa, and <gene>_outgroup2.fa."""


def write_two_genes(root: Path, genes: dict[str, dict[str, list[str]]] = TWO_GENES) -> TwoGeneDirs:
    """Write a gene set shaped like TWO_GENES under root in both layouts."""
    combined = root / "combined"
    separate = root / "separate"
    combined.mkdir()
    separate.mkdir()
    for gene, groups in genes.items():
        write_alignment_file(
            combined / f"{gene}.fa",
            {
                f"{gene}_{SPECIES[group]}_{i}": seq
                for group, seqs in groups.items()
                for i, seq in enumerate(seqs, start=1)
            },
        )
        for group, seqs in groups.items():
            write_alignment_file(
                separate / f"{gene}_{group}.fa",
                {f"{group}_{i}": seq for i, seq in enumerate(seqs, start=1)},
            )
    return TwoGeneDirs(combined=combined, separate=separate)
