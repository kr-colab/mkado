"""GFF3 annotation parser for extracting CDS regions."""

from __future__ import annotations

import gzip
import logging
import re
import zlib
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

from mkado.core.cds import CdsRegion
from mkado.errors import InputError

logger = logging.getLogger(__name__)


class GffFormatError(InputError):
    """A file given as GFF3 cannot be read as GFF3."""


def _parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3 attribute column (column 9)."""
    attrs: dict[str, str] = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, _, value = item.partition("=")
        attrs[key] = value
    return attrs


def _looks_like_gtf(attr_string: str) -> bool:
    """Detect GTF-style attributes (key "value") vs GFF3-style (key=value)."""
    # GTF uses: gene_id "ENSG00000..."; transcript_id "ENST00000...";
    # GFF3 uses: ID=gene1;Parent=tx1;Name=Foo
    return bool(re.search(r'\w+\s+"[^"]*"', attr_string))


def _numbered_lines(path: Path) -> Iterator[tuple[int, str]]:
    """Each line of a plain or gzipped GFF3 file with its line number, counted from 1.

    A file that cannot be read, is not text, or is not the gzip file its name says
    raises GffFormatError with the file name, not an error from deep inside the reader.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt") as f:
            yield from enumerate(f, 1)
    except FileNotFoundError:
        # A missing file keeps its standard error type, which callers test for.
        raise
    # gzip.BadGzipFile is an OSError; zlib.error comes from damaged compressed data.
    except (OSError, UnicodeDecodeError, EOFError, zlib.error) as e:
        raise GffFormatError(f"{path.name} cannot be read as GFF3: {e}") from e


def _check_transcript(cds_list: list[dict]) -> tuple[list[dict], int, int, str | None]:
    """A transcript's CDS lines sorted by start, the phase used, the duplicates removed, a reason.

    The CDS lines of one transcript are joined into one coding sequence and read in one
    frame. So they must lie on one contig and one strand, and must not overlap. The CDS
    first in coding order must give a valid phase, and the length after that phase must
    be a multiple of 3. The reason is None when the transcript passes; otherwise it says
    why it is skipped. Overlap with other transcripts or genes is not checked.

    A line that repeats another exactly is removed. Two copies of one interval that
    disagree on the phase are not one line: they overlap, so the transcript is skipped.
    """
    seen: set[tuple] = set()
    unique = []
    for cds in cds_list:
        key = (cds["chrom"], cds["start"], cds["end"], cds["strand"], cds["phase"])
        if key not in seen:
            seen.add(key)
            unique.append(cds)
    removed = len(cds_list) - len(unique)

    strands = {cds["strand"] for cds in unique}
    if len({cds["chrom"] for cds in unique}) > 1 or len(strands) > 1 or not strands <= {"+", "-"}:
        return [], 0, removed, "on more than one contig or strand"
    by_start = sorted(unique, key=lambda cds: cds["start"])
    if any(nxt["start"] < prev["end"] for prev, nxt in zip(by_start, by_start[1:])):
        return [], 0, removed, "with overlapping CDS"
    phase = (by_start[-1] if strands == {"-"} else by_start[0])["phase"]
    if phase is None:
        return [], 0, removed, "with an invalid phase"
    if (_cds_length(by_start) - phase) % 3:
        return [], 0, removed, "with CDS length not divisible by 3"
    return by_start, phase, removed, None


def _cds_length(cds_list: list[dict]) -> int:
    return sum(cds["end"] - cds["start"] for cds in cds_list)


def parse_gff3(
    path: str | Path,
    gene_ids: set[str] | None = None,
) -> list[CdsRegion]:
    """Parse a GFF3 file and extract CDS regions grouped by transcript.

    Selects the longest transcript per gene. CDS features are grouped by
    their Parent attribute (which should be mRNA/transcript IDs), then
    transcripts are grouped by their parent gene.

    Args:
        path: Path to GFF3 file.
        gene_ids: If provided, only return CDS regions for these gene IDs.

    Returns:
        List of CdsRegion objects, one per gene (longest transcript).
    """
    path = Path(path)

    # Collect CDS features keyed by parent transcript ID
    cds_by_transcript: dict[str, list[dict]] = defaultdict(list)
    # Map transcript ID -> gene ID (from mRNA/transcript features)
    transcript_to_gene: dict[str, str] = {}
    # Map gene ID -> gene name for display
    gene_names: dict[str, str] = {}
    # Track mRNA/transcript features
    transcript_info: dict[str, dict] = {}

    gtf_checked = False
    malformed_lines = 0
    for line_num, line in _numbered_lines(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 9:
            malformed_lines += 1
            logger.debug(
                "GFF3: skipping line %d (expected 9 tab-separated columns, got %d)",
                line_num,
                len(fields),
            )
            continue

        # Detect GTF format on first data line
        if not gtf_checked:
            gtf_checked = True
            if _looks_like_gtf(fields[8]):
                raise GffFormatError(
                    f"{path.name}: file appears to be GTF format, not GFF3 "
                    f"(line {line_num} has GTF-style attributes: "
                    f'key "value" instead of key=value). '
                    f"Convert to GFF3 or use a GTF-compatible parser."
                )

        # Parse coordinates with error handling
        try:
            start = int(fields[3]) - 1  # GFF3 is 1-based, convert to 0-based
            end = int(fields[4])  # GFF3 end is inclusive, but +1 for half-open
        except ValueError:
            malformed_lines += 1
            logger.debug(
                "GFF3: skipping line %d (non-integer coordinates: start=%r, end=%r)",
                line_num,
                fields[3],
                fields[4],
            )
            continue

        chrom = fields[0]
        feature_type = fields[2]
        strand = fields[6]
        phase_str = fields[7]
        attrs = _parse_attributes(fields[8])

        if feature_type == "gene":
            gid = attrs.get("ID", "")
            if gid:
                gene_names[gid] = attrs.get("Name", gid)

        elif feature_type in ("mRNA", "transcript"):
            tid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            if tid and parent:
                transcript_to_gene[tid] = parent
                transcript_info[tid] = {
                    "chrom": chrom,
                    "strand": strand,
                }

        elif feature_type == "CDS":
            parent = attrs.get("Parent", "")
            if not parent:
                logger.debug(
                    "GFF3: skipping CDS at line %d (no Parent attribute)",
                    line_num,
                )
                continue
            # Parent can be comma-separated in some GFF3 files
            for pid in parent.split(","):
                pid = pid.strip()
                # None marks a phase that is not valid; it matters only on the CDS whose
                # phase is used, which _check_transcript reads.
                phase = int(phase_str) if phase_str in ("0", "1", "2") else None
                cds_by_transcript[pid].append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": phase,
                    }
                )

    if malformed_lines > 0:
        logger.warning(
            "GFF3: %d malformed line(s) skipped (see debug log for details)",
            malformed_lines,
        )

    # Group transcripts by gene, select longest
    gene_transcripts: dict[str, list[str]] = defaultdict(list)
    for tid, gid in transcript_to_gene.items():
        if tid in cds_by_transcript:
            gene_transcripts[gid].append(tid)

    # Also handle CDS features that directly reference a gene (no mRNA level).
    # Per GFF3 spec NOTE 2, CDS may be parented directly to a gene.
    # We distinguish this (legitimate) from CDS parented to an unknown ID
    # (usually means mRNA features are missing or lack ID attributes).
    n_cds_to_gene = 0
    n_cds_to_unknown = 0
    for tid, cds_list in cds_by_transcript.items():
        if tid not in transcript_to_gene:
            if tid not in gene_transcripts:
                gene_transcripts[tid].append(tid)
                transcript_to_gene[tid] = tid
                if tid in gene_names:
                    n_cds_to_gene += 1
                else:
                    n_cds_to_unknown += 1

    if n_cds_to_unknown > 0:
        n_resolved = len(gene_transcripts) - n_cds_to_gene - n_cds_to_unknown
        if n_resolved > 0:
            logger.warning(
                "GFF3: %d CDS groups could not be linked to a gene via "
                "mRNA/transcript features and were treated as separate genes "
                "(only %d resolved normally). This usually means mRNA features "
                "are missing or lack an ID attribute — check your annotation file.",
                n_cds_to_unknown,
                n_resolved,
            )
        else:
            logger.warning(
                "GFF3: no mRNA/transcript features found — all %d CDS groups "
                "were treated as separate genes based on their Parent attribute. "
                "If your file has mRNA features, check that they have ID attributes.",
                n_cds_to_unknown,
            )

    results: list[CdsRegion] = []
    skipped_no_cds = 0
    skipped_no_valid = 0
    duplicates_removed = 0
    skipped_transcripts: dict[str, int] = defaultdict(int)

    for gene_id, tids in gene_transcripts.items():
        # Filter by gene IDs if specified
        display_name = gene_names.get(gene_id, gene_id)
        if gene_ids is not None:
            if gene_id not in gene_ids and display_name not in gene_ids:
                continue

        # Check each transcript first, so a gene with one bad transcript falls back to
        # its longest valid one.
        valid: dict[str, tuple[list[dict], int]] = {}
        for tid in tids:
            by_start, phase, removed, reason = _check_transcript(cds_by_transcript[tid])
            duplicates_removed += removed
            if reason is None:
                valid[tid] = (by_start, phase)
            else:
                skipped_transcripts[reason] += 1
                logger.debug("GFF3: skipping transcript %s of %s (CDS %s)", tid, gene_id, reason)
        if not valid:
            skipped_no_valid += 1
            logger.debug("GFF3: skipping %s (no valid transcript)", gene_id)
            continue

        # Select the longest valid transcript; the first one wins a tie.
        best_tid, (by_start, phase) = max(valid.items(), key=lambda item: _cds_length(item[1][0]))
        if _cds_length(by_start) <= 0:
            skipped_no_cds += 1
            logger.debug("GFF3: skipping %s (no transcript with CDS features)", gene_id)
            continue

        region = CdsRegion(
            gene_id=gene_id,
            transcript_id=best_tid,
            chrom=by_start[0]["chrom"],
            exons=[(c["start"], c["end"]) for c in by_start],
            strand=by_start[0]["strand"],
            phase=phase,
        )
        results.append(region)

    if duplicates_removed > 0:
        logger.warning("GFF3: %d duplicate CDS line(s) removed", duplicates_removed)
    if skipped_transcripts:
        logger.warning(
            "GFF3: %d transcript(s) skipped — %s",
            sum(skipped_transcripts.values()),
            "; ".join(f"{n} {reason}" for reason, n in skipped_transcripts.items()),
        )

    n_skipped = skipped_no_cds + skipped_no_valid
    if n_skipped > 0:
        parts = []
        if skipped_no_cds > 0:
            parts.append(f"{skipped_no_cds} with no CDS features")
        if skipped_no_valid > 0:
            parts.append(f"{skipped_no_valid} with no valid transcript")
        logger.warning(
            "GFF3: %d genes skipped — %s",
            n_skipped,
            "; ".join(parts),
        )

    return results
