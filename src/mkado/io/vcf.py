"""VCF polymorphism and divergence extraction for MK tests."""

from __future__ import annotations

import logging
import os
import re
import tempfile
import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mkado.analysis.asymptotic import PolymorphismData
from mkado.core.cds import CdsRegion, _COMPLEMENT
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.data.genetic_codes import CLEAN_CODONS

if TYPE_CHECKING:
    from collections.abc import Iterator

logger = logging.getLogger(__name__)

# An htslib log line: a level letter, the function that spoke, then the message
# to the end of the line. Matched wherever it starts, since another writer may
# have left the line unterminated just before it.
_HTSLIB_LINE = re.compile(rb"\[[EWIDT]::\w+\] [^\n]*\n?")

_capturing = False


@contextmanager
def _htslib_stderr_to_log() -> Iterator[None]:
    """Route what htslib writes to file descriptor 2 through the logger at DEBUG.

    htslib writes straight to fd 2, bypassing Python's warning system, both
    while a file is opened and while a region query parses a record, so a
    caller that wants the query-time messages holds this around the whole
    extraction. Inside the block fd 2 is a temporary file, which has no size
    ceiling where a pipe would block the writer once full. The redirect is
    process-wide, so a progress bar redraw from another thread lands in the
    file too: on exit, htslib's lines go to the logger and every other byte is
    written back to fd 2, delayed by the block. A block entered inside another
    is a no-op, so nested captures drain once, at the outer level.
    """
    global _capturing
    if _capturing:
        yield
        return
    _capturing = True
    try:
        with tempfile.TemporaryFile() as capture:
            orig_fd = os.dup(2)
            os.dup2(capture.fileno(), 2)
            try:
                yield
            finally:
                os.dup2(orig_fd, 2)
                os.close(orig_fd)
                capture.seek(0)
                _drain(capture.read())
    finally:
        _capturing = False


def _drain(data: bytes) -> None:
    """Log htslib's lines wherever they start and write everything else back to fd 2."""
    end = 0
    for match in _HTSLIB_LINE.finditer(data):
        if match.start() > end:
            os.write(2, data[end : match.start()])
        logger.debug("htslib: %s", match.group().decode(errors="replace").rstrip("\r\n"))
        end = match.end()
    if end < len(data):
        os.write(2, data[end:])


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


def genotyped_samples(path: str | Path) -> list[str]:
    """Sample names of a VCF that carries genotypes, or an empty list if it carries none.

    A VCF carries no genotypes when it has no sample columns, or when its header
    declares no ``GT`` FORMAT field. No result can come from such a file: an
    ingroup gives no polymorphism, and an outgroup reads as the reference base
    everywhere.
    """
    vcf = _open_vcf(path)
    try:
        has_gt = any(
            header.type == "FORMAT" and header.info().get("ID") == "GT"
            for header in vcf.header_iter()
        )
        return list(vcf.samples) if has_gt else []
    finally:
        vcf.close()


@dataclass
class _SnpInfo:
    """A biallelic ingroup SNP with its called allele counts. ``n_alt`` is at least 1."""

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
    skipped_ambiguous: int = 0
    """SNPs and fixed differences dropped because a codon has a base other than A, C, G or T."""
    skipped_filtered: int = 0
    """Ingroup records with a FILTER value other than PASS or '.'."""
    skipped_ref_mismatch: int = 0
    """Ingroup records whose REF is not the reference FASTA base, and ingroup SNPs in a
    codon that holds such a position in either VCF."""
    outgroup_skipped: int = 0
    """Outgroup records dropped for any reason, including a record with no genotype."""

    def count(self, reason: _Skip) -> None:
        """Add one ingroup record skipped for ``reason``."""
        match reason:
            case _Skip.INDEL:
                self.skipped_indels += 1
            case _Skip.MULTIALLELIC:
                self.skipped_multiallelic += 1
            case _Skip.FILTERED:
                self.skipped_filtered += 1
            case _Skip.REF_MISMATCH:
                self.skipped_ref_mismatch += 1
            case _Skip.INVARIANT:
                # No ALT allele: an invariant site, so no data is lost.
                pass


class _Skip(Enum):
    """Why a VCF record is not a usable biallelic SNV."""

    INVARIANT = auto()
    INDEL = auto()
    MULTIALLELIC = auto()
    FILTERED = auto()
    REF_MISMATCH = auto()


def _skip_reason(
    variant: object, ref_fetch: callable, chrom: str, include_filtered: bool
) -> _Skip | None:
    """Why a record cannot be used as a biallelic SNV, or None if it can."""
    ref = variant.REF.upper()
    # Checked first: a mismatch says the FASTA base is wrong, whatever else the record is.
    if len(ref) == 1 and ref != ref_fetch(chrom, variant.POS - 1):
        return _Skip.REF_MISMATCH
    if not variant.ALT:
        return _Skip.INVARIANT
    if variant.is_indel:
        return _Skip.INDEL
    if len(variant.ALT) != 1:
        return _Skip.MULTIALLELIC
    if not _is_snv(ref, variant.ALT[0]):
        return _Skip.INDEL
    # cyvcf2 gives None for both PASS and '.'.
    if variant.FILTER is not None and not include_filtered:
        return _Skip.FILTERED
    return None


def _is_snv(ref: str, alt: str) -> bool:
    """Whether a biallelic record is a single-base substitution.

    A ``*`` ALT is one character, but it marks a deletion that spans the site,
    so it is not a base.
    """
    return len(ref) == 1 and len(alt) == 1 and alt != "*"


def _genotype_alleles(variant: object) -> np.ndarray | None:
    """Allele indices of each sample's call (samples x copies), or None with no samples.

    A missing allele is -1, and a copy past a sample's own ploidy, such as a
    haploid call in a diploid record, is -2. cyvcf2's ``gt_types`` is not used
    because it types ``0/.`` as HOM_REF and ``./1`` as HET, and so adds copies
    the call does not have.
    """
    genotype = variant.genotype
    if genotype is None:
        return None
    # The last column is the phase flag, not an allele.
    return genotype.array()[:, :-1]


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


def _query_ingroup_snps(
    vcf: object,
    cds: CdsRegion,
    stats: GeneStats,
    ref_fetch: callable,
    include_filtered: bool,
) -> tuple[list[_SnpInfo], set[int]]:
    """Query ingroup VCF for biallelic SNPs overlapping a CDS region.

    Returns the SNPs in genomic coordinates with their allele counts, and the
    positions whose REF is not the FASTA base.
    """
    snps: list[_SnpInfo] = []
    mismatched: set[int] = set()

    # Build query regions from exons
    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"  # cyvcf2 uses 1-based
        for variant in _query_region(vcf, region):
            pos_0 = variant.POS - 1  # convert to 0-based

            # Must be in our CDS
            if not cds.contains_position(pos_0):
                continue

            reason = _skip_reason(variant, ref_fetch, cds.chrom, include_filtered)
            if reason is _Skip.REF_MISMATCH:
                mismatched.add(pos_0)
            if reason is not None:
                stats.count(reason)
                continue
            ref = variant.REF.upper()
            alt = variant.ALT[0].upper()

            # Count each called allele. The site is biallelic, so a called allele is
            # 0 (REF) or 1 (ALT), and a missing one adds nothing.
            alleles = _genotype_alleles(variant)
            n_ref = n_alt = 0
            if alleles is not None:
                n_ref = int(np.count_nonzero(alleles == 0))
                n_alt = int(np.count_nonzero(alleles == 1))

            total = n_ref + n_alt
            if total == 0:
                stats.skipped_missing += 1
                continue
            # No called ALT allele: the same as no record, so the ingroup carries REF.
            if n_alt == 0:
                continue

            snps.append(_SnpInfo(pos=pos_0, ref=ref, alt=alt, n_ref=n_ref, n_alt=n_alt))

    return snps, mismatched


def _query_outgroup_genotype(
    vcf: object,
    cds: CdsRegion,
    stats: GeneStats,
    ref_fetch: callable,
    include_filtered: bool,
) -> tuple[dict[int, str], set[int]]:
    """Query outgroup VCF for genotypes at CDS positions.

    Returns a dict mapping 0-based genomic position -> outgroup allele, and the
    positions whose REF is not the FASTA base. The dict only includes positions
    where every allele of the outgroup's call is called and is ALT. A
    heterozygous or partly missing call cannot resolve a single outgroup allele,
    so it is left out and treated like a missing record.
    """
    outgroup_alleles: dict[int, str] = {}
    mismatched: set[int] = set()

    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"
        for variant in _query_region(vcf, region):
            pos_0 = variant.POS - 1

            if not cds.contains_position(pos_0):
                continue

            reason = _skip_reason(variant, ref_fetch, cds.chrom, include_filtered)
            if reason is _Skip.REF_MISMATCH:
                mismatched.add(pos_0)
            if reason is not None:
                if reason is not _Skip.INVARIANT:
                    stats.outgroup_skipped += 1
                continue
            alleles = _genotype_alleles(variant)
            if alleles is None:
                stats.outgroup_skipped += 1
                continue
            alt = variant.ALT[0].upper()

            # The outgroup is the first sample. Drop the padding past its ploidy.
            call = alleles[0][alleles[0] != -2]
            if call.size and (call == 1).all():
                outgroup_alleles[pos_0] = alt

    return outgroup_alleles, mismatched


def extract_gene_data(
    vcf_path: str | Path,
    outgroup_vcf_path: str | Path | None,
    cds: CdsRegion,
    ref_fasta_path: str | Path,
    genetic_code: GeneticCode | None = None,
    min_frequency: float = 0.0,
    no_singletons: bool = False,
    *,
    include_filtered: bool = False,
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
        include_filtered: If True, also use records whose FILTER is not PASS or '.'.
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

    # A CDS that runs past its contig means the annotation does not match the reference.
    contig_length = ref_fasta.get_reference_length(cds.chrom)
    cds_end = cds.exons[-1][1]
    if cds_end > contig_length:
        if owns_ref_fasta:
            ref_fasta.close()
        raise ValueError(f"CDS ends at {cds_end}, past the end of {cds.chrom} ({contig_length} bp)")
    ref_fetch = _ref_base_fetcher(ref_fasta)

    # Open each VCF that was not passed in, and close it after its query.
    vcf = ingroup_vcf or _open_vcf(vcf_path)
    ingroup_snps, mismatched = _query_ingroup_snps(vcf, cds, stats, ref_fetch, include_filtered)
    if ingroup_vcf is None:
        vcf.close()

    outgroup_alleles: dict[int, str] = {}
    if outgroup_vcf is not None or outgroup_vcf_path is not None:
        vcf = outgroup_vcf or _open_vcf(outgroup_vcf_path)
        outgroup_alleles, outgroup_mismatched = _query_outgroup_genotype(
            vcf, cds, stats, ref_fetch, include_filtered
        )
        mismatched |= outgroup_mismatched
        if outgroup_vcf is None:
            vcf.close()

    # The FASTA base at a REF mismatch cannot be trusted, so no record in a codon that
    # holds one is used, in either group.
    mismatch_codons = {cds.genomic_pos_to_codon_index(pos) for pos in mismatched}
    if mismatch_codons:

        def usable(pos: int) -> bool:
            return cds.genomic_pos_to_codon_index(pos) not in mismatch_codons

        kept_snps = [s for s in ingroup_snps if usable(s.pos)]
        stats.skipped_ref_mismatch += len(ingroup_snps) - len(kept_snps)
        ingroup_snps = kept_snps
        kept_alleles = {pos: alt for pos, alt in outgroup_alleles.items() if usable(pos)}
        stats.outgroup_skipped += len(outgroup_alleles) - len(kept_alleles)
        outgroup_alleles = kept_alleles

    # === Polymorphism extraction ===
    polymorphisms: list[tuple[float, str]] = []

    for snp in ingroup_snps:
        # Skip sites fixed for ALT in ingroup (these are potential divergence, not polymorphism)
        if snp.is_fixed_alt:
            continue

        codon_idx = cds.genomic_pos_to_codon_index(snp.pos)
        if codon_idx is None:
            continue

        # Reconstruct codons
        ref_codon, alt_codon = _reconstruct_codon_with_sub(
            cds, codon_idx, ref_fetch, snp.pos, snp.alt
        )

        # An ambiguous reference base (N, an IUPAC code) leaves the change unknown.
        if not {ref_codon, alt_codon} <= CLEAN_CODONS:
            stats.skipped_ambiguous += 1
            continue

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
            # Checked before get_path, which has no path for such a codon and would
            # otherwise report the pair as stop-blocked.
            if not {in_codon, out_codon} <= CLEAN_CODONS:
                stats.skipped_ambiguous += 1
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
