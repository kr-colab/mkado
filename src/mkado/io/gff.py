"""GFF3 annotation parser for extracting CDS regions."""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

from mkado.core.cds import CdsRegion


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

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3]) - 1  # GFF3 is 1-based, convert to 0-based
            end = int(fields[4])  # GFF3 end is inclusive, but +1 for half-open
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
                    continue
                # Parent can be comma-separated in some GFF3 files
                for pid in parent.split(","):
                    pid = pid.strip()
                    phase = int(phase_str) if phase_str in ("0", "1", "2") else 0
                    cds_by_transcript[pid].append(
                        {
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "phase": phase,
                        }
                    )

    # Group transcripts by gene, select longest
    gene_transcripts: dict[str, list[str]] = defaultdict(list)
    for tid, gid in transcript_to_gene.items():
        if tid in cds_by_transcript:
            gene_transcripts[gid].append(tid)

    # Also handle CDS features that directly reference a gene (no mRNA level)
    for tid, cds_list in cds_by_transcript.items():
        if tid not in transcript_to_gene:
            # This CDS parent is a gene directly, or an unknown parent
            # Treat tid as both gene and transcript
            if tid not in gene_transcripts:
                gene_transcripts[tid].append(tid)
                transcript_to_gene[tid] = tid

    results: list[CdsRegion] = []

    for gene_id, tids in gene_transcripts.items():
        # Filter by gene IDs if specified
        display_name = gene_names.get(gene_id, gene_id)
        if gene_ids is not None:
            if gene_id not in gene_ids and display_name not in gene_ids:
                continue

        # Select longest transcript
        best_tid = None
        best_length = 0
        for tid in tids:
            cds_list = cds_by_transcript[tid]
            length = sum(c["end"] - c["start"] for c in cds_list)
            if length > best_length:
                best_length = length
                best_tid = tid

        if best_tid is None:
            continue

        cds_list = cds_by_transcript[best_tid]
        if not cds_list:
            continue

        exons = [(c["start"], c["end"]) for c in cds_list]
        exons.sort(key=lambda e: e[0])

        chrom = cds_list[0]["chrom"]
        strand = cds_list[0]["strand"]

        # Phase of the first CDS in coding order
        sorted_cds = sorted(cds_list, key=lambda c: c["start"])
        if strand == "-":
            first_cds_phase = sorted_cds[-1]["phase"]
        else:
            first_cds_phase = sorted_cds[0]["phase"]

        region = CdsRegion(
            gene_id=gene_id,
            transcript_id=best_tid,
            chrom=chrom,
            exons=exons,
            strand=strand,
            phase=first_cds_phase,
        )

        # Skip genes where CDS length is not divisible by 3
        if not region.is_valid():
            continue

        results.append(region)

    return results
