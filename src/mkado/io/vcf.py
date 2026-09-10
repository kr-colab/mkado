"""VCF polymorphism and divergence extraction for MK tests."""

from __future__ import annotations

import logging
import os
import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from mkado.analysis.asymptotic import PolymorphismData
from mkado.core.cds import CdsRegion, _COMPLEMENT
from mkado.core.codons import DEFAULT_CODE, GeneticCode

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)


@contextmanager
def _htslib_stderr_to_log() -> Iterator[None]:
    """Route what htslib writes to file descriptor 2 through the logger at DEBUG.

    htslib writes warnings (e.g., "FORMAT 'GT' not defined in header") and index
    errors directly to fd 2, bypassing Python's warning system. Inside this
    block fd 2 is a pipe, drained into the logger afterwards.
    """
    r_fd, w_fd = os.pipe()
    orig_fd = os.dup(2)
    os.dup2(w_fd, 2)
    os.close(w_fd)
    try:
        yield
    finally:
        # Restore fd 2 before draining: read() only returns once the pipe's last writer is gone.
        os.dup2(orig_fd, 2)
        os.close(orig_fd)
        with os.fdopen(r_fd) as f:
            for line in f.read().strip().splitlines():
                logger.debug("htslib: %s", line)


def _open_vcf(path: str | Path) -> object:
    """Open a cyvcf2.VCF, capturing htslib stderr warnings via the Python logger."""
    import cyvcf2

    with _htslib_stderr_to_log():
        return cyvcf2.VCF(str(path))


class VcfQueryError(RuntimeError):
    """A VCF could not be opened or could not answer a region query."""


def _query_region(vcf: object, region: str) -> list:
    """Return the records in a region, raising when the file cannot be queried.

    A contig the file does not carry returns nothing with a Python warning,
    which is a legitimate empty result, so that warning is silenced. The query
    itself raises only when the index is unusable, and an unanswered query must
    not pass as "no variants here".
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return list(vcf(region))
        except Exception as exc:
            raise VcfQueryError(f"region query {region} failed: {exc}") from exc


def probe_vcf(path: str | Path) -> list[str]:
    """Open a VCF, prove it answers a region query, and return its header contigs.

    Raises:
        VcfQueryError: If the file cannot be opened, or has no usable index.
    """
    try:
        vcf = _open_vcf(path)
    except Exception as exc:
        raise VcfQueryError(f"{path} cannot be opened: {exc}") from exc
    try:
        with _htslib_stderr_to_log():
            try:
                contigs = list(vcf.seqnames)
                # Any contig name proves the index: through an index a name the
                # file lacks returns no records, and without one the query raises.
                _query_region(vcf, f"{contigs[0] if contigs else 'probe'}:1-1")
            except Exception as exc:
                raise VcfQueryError(
                    f"{path} has no usable index. Bgzip the file and index it with "
                    f"tabix or bcftools index. ({exc})"
                ) from exc
    finally:
        vcf.close()
    return contigs


@dataclass
class _SnpInfo:
    """A biallelic ingroup SNP with its diploid allele counts."""

    pos: int  # 0-based genomic position
    ref: str
    alt: str
    n_ref: int
    n_alt: int

    @property
    def total(self) -> int:
        return self.n_ref + self.n_alt

    @property
    def alt_freq(self) -> float:
        return self.n_alt / self.total

    @property
    def is_fixed_alt(self) -> bool:
        return self.n_ref == 0


@dataclass
class GeneStats:
    """Bookkeeping stats for a single gene's VCF extraction."""

    skipped_indels: int = 0
    skipped_multiallelic: int = 0
    skipped_missing: int = 0
    stop_blocked_pairs: int = 0
    """Fixed differences with no stop-free mutational ordering, dropped from Dn/Ds."""


def _ref_base_fetcher(fasta_file: object) -> callable:
    """Create a reference base fetcher from a pysam.FastaFile."""

    def fetch(chrom: str, pos: int) -> str:
        return fasta_file.fetch(chrom, pos, pos + 1).upper()

    return fetch


def _complement_base(base: str) -> str:
    """Complement a single base."""
    return base.translate(_COMPLEMENT)


def _reconstruct_codon_with_sub(
    cds: CdsRegion,
    codon_index: int,
    ref_fetch: callable,
    sub_pos: int,
    sub_base: str,
) -> tuple[str, str]:
    """Reconstruct ref and alt codons given a SNP.

    Args:
        cds: CDS region.
        codon_index: Codon index in CDS.
        ref_fetch: Callable to get reference base at (chrom, pos).
        sub_pos: Genomic position of the substitution.
        sub_base: The alternative base (in forward-strand terms).

    Returns:
        (ref_codon, alt_codon) in coding-strand orientation.
    """
    p1, p2, p3 = cds.codon_positions(codon_index)
    bases = [ref_fetch(cds.chrom, p1), ref_fetch(cds.chrom, p2), ref_fetch(cds.chrom, p3)]

    # Determine which position in the codon the SNP falls
    offset = cds.genomic_pos_to_codon_offset(sub_pos)
    alt_bases = list(bases)
    alt_bases[offset] = sub_base

    ref_codon = "".join(bases)
    alt_codon = "".join(alt_bases)

    if cds.strand == "-":
        ref_codon = ref_codon.translate(_COMPLEMENT)
        alt_codon = alt_codon.translate(_COMPLEMENT)

    return ref_codon.upper(), alt_codon.upper()


def _query_ingroup_snps_with_handle(
    vcf: object,
    cds: CdsRegion,
    stats: GeneStats,
) -> list[_SnpInfo]:
    """Query ingroup VCF for biallelic SNPs overlapping a CDS region.

    Uses a pre-opened cyvcf2.VCF handle. Returns SNPs in genomic coordinates
    with allele frequencies.
    """
    snps: list[_SnpInfo] = []

    # Build query regions from exons
    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"  # cyvcf2 uses 1-based
        for variant in _query_region(vcf, region):
            pos_0 = variant.POS - 1  # convert to 0-based

            # Must be in our CDS
            if not cds.contains_position(pos_0):
                continue

            # Skip indels
            if variant.is_indel:
                stats.skipped_indels += 1
                continue

            # Skip multi-allelic
            if len(variant.ALT) != 1:
                stats.skipped_multiallelic += 1
                continue

            ref = variant.REF.upper()
            alt = variant.ALT[0].upper()

            # Must be SNP
            if len(ref) != 1 or len(alt) != 1:
                stats.skipped_indels += 1
                continue

            # Get allele frequency from genotypes
            gt_types = variant.gt_types  # 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
            n_ref = 0
            n_alt = 0
            for gt in gt_types:
                if gt == 0:  # HOM_REF
                    n_ref += 2
                elif gt == 1:  # HET
                    n_ref += 1
                    n_alt += 1
                elif gt == 3:  # HOM_ALT
                    n_alt += 2
                # gt == 2 is missing, skip

            total = n_ref + n_alt
            if total == 0:
                stats.skipped_missing += 1
                continue

            snps.append(_SnpInfo(pos=pos_0, ref=ref, alt=alt, n_ref=n_ref, n_alt=n_alt))

    return snps


def _query_ingroup_snps(
    vcf_path: str | Path,
    cds: CdsRegion,
    stats: GeneStats,
) -> list[_SnpInfo]:
    """Query ingroup VCF for biallelic SNPs overlapping a CDS region."""
    vcf = _open_vcf(vcf_path)
    snps = _query_ingroup_snps_with_handle(vcf, cds, stats)
    vcf.close()
    return snps


def _query_outgroup_genotype_with_handle(
    vcf: object,
    cds: CdsRegion,
) -> dict[int, str]:
    """Query outgroup VCF for genotypes at CDS positions.

    Uses a pre-opened cyvcf2.VCF handle. Returns dict mapping 0-based genomic
    position -> outgroup allele. Only includes positions where the outgroup is
    homozygous for an allele that differs from reference; a heterozygous call
    cannot resolve a single outgroup allele, so it is left out and treated like
    a missing record.
    """
    outgroup_alleles: dict[int, str] = {}

    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"
        for variant in _query_region(vcf, region):
            pos_0 = variant.POS - 1

            if not cds.contains_position(pos_0):
                continue

            # Skip indels and multi-allelic
            if variant.is_indel or len(variant.ALT) != 1:
                continue

            ref = variant.REF.upper()
            alt = variant.ALT[0].upper()

            if len(ref) != 1 or len(alt) != 1:
                continue

            # For single-sample outgroup, check if it carries ALT
            gt_types = variant.gt_types
            if len(gt_types) == 0:
                continue

            gt = gt_types[0]
            if gt == 3:  # HOM_ALT
                outgroup_alleles[pos_0] = alt

    return outgroup_alleles


def _query_outgroup_genotype(
    vcf_path: str | Path,
    cds: CdsRegion,
) -> dict[int, str]:
    """Query outgroup VCF for genotypes at CDS positions."""
    vcf = _open_vcf(vcf_path)
    outgroup_alleles = _query_outgroup_genotype_with_handle(vcf, cds)
    vcf.close()
    return outgroup_alleles


def extract_gene_data(
    vcf_path: str | Path,
    outgroup_vcf_path: str | Path | None,
    cds: CdsRegion,
    ref_fasta_path: str | Path,
    genetic_code: GeneticCode | None = None,
    min_frequency: float = 0.0,
    no_singletons: bool = False,
    *,
    ingroup_vcf: object | None = None,
    outgroup_vcf: object | None = None,
    ref_fasta: object | None = None,
) -> tuple[PolymorphismData, GeneStats]:
    """Extract polymorphism and divergence data for a single gene from VCF.

    Args:
        vcf_path: Path to ingroup VCF file.
        outgroup_vcf_path: Path to outgroup VCF file, or None.
        cds: CDS region defining the gene.
        ref_fasta_path: Path to indexed reference FASTA.
        genetic_code: Genetic code for classification (standard if None).
        min_frequency: Minimum derived allele frequency threshold.
        no_singletons: If True, exclude sites where the derived allele is seen once.
        ingroup_vcf: Pre-opened cyvcf2.VCF handle for ingroup (optional).
        outgroup_vcf: Pre-opened cyvcf2.VCF handle for outgroup (optional).
        ref_fasta: Pre-opened pysam.FastaFile handle (optional).

    Returns:
        Tuple of (PolymorphismData, GeneStats).
    """
    code = genetic_code or DEFAULT_CODE
    stats = GeneStats()

    # Use pre-opened handles if provided, otherwise open new ones
    owns_ref_fasta = ref_fasta is None
    if owns_ref_fasta:
        import pysam

        ref_fasta = pysam.FastaFile(str(ref_fasta_path))
    ref_fetch = _ref_base_fetcher(ref_fasta)

    # Get ingroup SNPs
    if ingroup_vcf is not None:
        ingroup_snps = _query_ingroup_snps_with_handle(ingroup_vcf, cds, stats)
    else:
        ingroup_snps = _query_ingroup_snps(vcf_path, cds, stats)

    # Get outgroup genotypes
    outgroup_alleles: dict[int, str] = {}
    if outgroup_vcf is not None:
        outgroup_alleles = _query_outgroup_genotype_with_handle(outgroup_vcf, cds)
    elif outgroup_vcf_path is not None:
        outgroup_alleles = _query_outgroup_genotype(outgroup_vcf_path, cds)

    # === Polymorphism extraction ===
    polymorphisms: list[tuple[float, str]] = []

    for snp in ingroup_snps:
        # Skip sites fixed for ALT in ingroup (these are potential divergence, not polymorphism)
        if snp.is_fixed_alt:
            continue

        # Skip sites fixed for REF (freq == 0)
        if snp.alt_freq <= 0.0 or snp.alt_freq >= 1.0:
            continue

        codon_idx = cds.genomic_pos_to_codon_index(snp.pos)
        if codon_idx is None:
            continue

        # Reconstruct codons
        ref_codon, alt_codon = _reconstruct_codon_with_sub(
            cds, codon_idx, ref_fetch, snp.pos, snp.alt
        )

        # Skip stop codons in reference
        if code.translate(ref_codon) == "*":
            continue

        # Classify the change
        is_syn = code.is_synonymous_change(ref_codon, alt_codon)
        if is_syn is None:
            continue

        # If the outgroup carries ALT, REF is the derived allele.
        derived_count = snp.n_alt
        if snp.pos in outgroup_alleles and outgroup_alleles[snp.pos] == snp.alt:
            derived_count = snp.n_ref
        derived_freq = derived_count / snp.total

        # A singleton is one copy of the derived allele, whatever the number of
        # called samples, so it is filtered by count rather than by frequency.
        if no_singletons and derived_count <= 1:
            continue
        if derived_freq < min_frequency:
            continue
        if derived_freq <= 0.0 or derived_freq >= 1.0:
            continue

        poly_type = "S" if is_syn else "N"
        polymorphisms.append((derived_freq, poly_type))

    # === Divergence extraction ===
    dn = 0
    ds = 0

    if outgroup_vcf_path is not None or outgroup_vcf is not None:
        # The ingroup codon carries ALT where the ingroup is fixed for it. A codon needs
        # one clean codon per group to be a fixed difference, so any polymorphic position,
        # including a third-allele outgroup base there, voids the whole codon, not just
        # that one base. A position with no outgroup record is the reference.
        fixed_alt = {s.pos: s.alt for s in ingroup_snps if s.is_fixed_alt}
        polymorphic = {s.pos for s in ingroup_snps if not s.is_fixed_alt}
        poly_codons = {cds.genomic_pos_to_codon_index(pos) for pos in polymorphic}

        def in_fetch(chrom: str, pos: int) -> str:
            return fixed_alt.get(pos) or ref_fetch(chrom, pos)

        def out_fetch(chrom: str, pos: int) -> str:
            return outgroup_alleles.get(pos) or ref_fetch(chrom, pos)

        # Both allele maps were filtered by cds.contains_position, so every lookup hits.
        candidates = {
            cds.genomic_pos_to_codon_index(pos)
            for pos in fixed_alt.keys() | outgroup_alleles.keys()
        } - poly_codons
        for codon_idx in sorted(candidates):
            in_codon = cds.extract_codon(codon_idx, in_fetch)
            out_codon = cds.extract_codon(codon_idx, out_fetch)
            if in_codon == out_codon:
                continue
            if code.translate(in_codon) == "*" or code.translate(out_codon) == "*":
                continue
            # get_path handles codons that differ at more than one position.
            path = code.get_path(in_codon, out_codon)
            if not path:
                stats.stop_blocked_pairs += 1
                continue
            for change_type, _position in path:
                if change_type == "R":
                    dn += 1
                elif change_type == "S":
                    ds += 1

    if owns_ref_fasta:
        ref_fasta.close()

    return PolymorphismData(
        polymorphisms=polymorphisms,
        dn=dn,
        ds=ds,
        gene_id=cds.gene_id,
    ), stats
