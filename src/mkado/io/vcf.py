"""VCF polymorphism and divergence extraction for MK tests."""

from __future__ import annotations

import logging
import os
import re
import tempfile
import warnings
from collections import Counter, defaultdict
from itertools import permutations, product
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mkado.analysis.asymptotic import PolymorphismData
from mkado.core.cds import CdsRegion, _COMPLEMENT
from mkado.core.alignment import passes_frequency_filter
from mkado.core.codons import DEFAULT_CODE, GeneticCode
from mkado.core.sequences import ancestral_base
from mkado.data.genetic_codes import CLEAN_CODONS
from mkado.errors import InputError

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


def probe_ref(path: str | Path) -> None:
    """Open a reference FASTA once, so that a missing index is built before any worker opens it.

    pysam builds a missing .fai (and .gzi) index when it opens a FASTA. Workers that each
    open an unindexed reference race to write the same index, and the genes they read
    from a part-written index fail.

    Raises:
        InputError: If the file cannot be opened as an indexed FASTA.
    """
    import pysam

    try:
        # htslib writes its own error line to fd 2; the raised message says the same.
        with _htslib_stderr_to_log():
            pysam.FastaFile(str(path)).close()
    except OSError as exc:
        raise InputError(f"{path} cannot be opened as an indexed FASTA: {exc}") from exc


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
    skipped_stop: int = 0
    """SNPs and fixed differences dropped because a codon is a stop codon."""
    skipped_ambiguous: int = 0
    """SNPs and fixed differences dropped because a codon has a base other than A, C, G or T."""
    skipped_filtered: int = 0
    """Ingroup records with a FILTER value other than PASS or '.'."""
    skipped_ref_mismatch: int = 0
    """Ingroup records whose REF is not the reference FASTA base, and ingroup SNPs in a
    codon that holds such a position in either VCF."""
    outgroup_skipped: int = 0
    """Outgroup records dropped for any reason, including a record with no genotype."""
    swapped_ref_alt: int = 0
    """Records in either VCF used with REF and ALT swapped: the VCF REF was not the FASTA
    base, but the VCF ALT was."""
    skipped_unphased: int = 0
    """Codons whose fixed difference depends on the unknown phase of two heterozygous
    outgroup calls."""

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


class _CodonSkip(Enum):
    """Why a codon pair adds no fixed difference, where that is counted."""

    AMBIGUOUS = auto()
    STOP = auto()
    STOP_BLOCKED = auto()


class _SnpDrop(Enum):
    """Why one outgroup phase drops a polymorphism."""

    MISSING = auto()
    """No outgroup codon is clean."""
    UNPOLARIZED = auto()
    """No base is shared with the outgroup, and a derived frequency is needed."""
    FILTERED = auto()
    """The frequency filters remove it."""


def _in_full_codon(cds: CdsRegion, pos: int) -> bool:
    """Whether ``pos`` is in a full codon of ``cds``; bases after the last one are in none."""
    codon_idx = cds.genomic_pos_to_codon_index(pos)
    return codon_idx is not None and codon_idx < cds.num_codons()


class _Skip(Enum):
    """Why a VCF record is not a usable biallelic SNV."""

    INVARIANT = auto()
    INDEL = auto()
    MULTIALLELIC = auto()
    FILTERED = auto()
    REF_MISMATCH = auto()


def _skip_reason(variant: object, fasta_base: str, include_filtered: bool) -> _Skip | None:
    """Why a record cannot be used as a biallelic SNV, or None if it can.

    A biallelic SNV whose REF is not ``fasta_base`` but whose ALT is can be used:
    ``_snv_alleles`` reads it with the two alleles swapped.
    """
    ref = variant.REF.upper()
    # Checked first: unless the ALT is the FASTA base (the alleles are only swapped), a
    # mismatch says the FASTA base is wrong, whatever else the record is.
    alt_is_fasta = len(variant.ALT) == 1 and variant.ALT[0].upper() == fasta_base
    if len(ref) == 1 and ref != fasta_base and not alt_is_fasta:
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


def _snv_alleles(variant: object, fasta_base: str) -> tuple[np.ndarray | None, str, str, bool]:
    """Allele indices, REF, ALT and swap flag of a record that ``_skip_reason`` passed.

    When the VCF REF is not the FASTA base, the ALT is, so the two are swapped back:
    REF and ALT trade places and called indices 0 and 1 are exchanged. Index 0 is
    then always the FASTA base, and index 1 the other base.
    """
    ref = variant.REF.upper()
    alt = variant.ALT[0].upper()
    alleles = _genotype_alleles(variant)
    swapped = ref != fasta_base
    if swapped:
        ref, alt = alt, ref
        if alleles is not None:
            alleles = np.where(alleles >= 0, 1 - alleles, alleles)
    return alleles, ref, alt, swapped


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
) -> tuple[list[_SnpInfo], set[int], set[int]]:
    """Query ingroup VCF for biallelic SNPs overlapping a CDS region.

    Returns the SNPs in genomic coordinates with their allele counts, the
    positions whose REF is not the FASTA base, and the positions of the SNPs
    read with REF and ALT swapped.
    """
    snps: list[_SnpInfo] = []
    mismatched: set[int] = set()
    swapped_at: set[int] = set()

    # Build query regions from exons
    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"  # cyvcf2 uses 1-based
        for variant in _query_region(vcf, region):
            pos_0 = variant.POS - 1  # convert to 0-based

            # Must be in our CDS
            if not _in_full_codon(cds, pos_0):
                continue

            fasta_base = ref_fetch(cds.chrom, pos_0)
            reason = _skip_reason(variant, fasta_base, include_filtered)
            if reason is _Skip.REF_MISMATCH:
                mismatched.add(pos_0)
            if reason is not None:
                stats.count(reason)
                continue
            alleles, ref, alt, swapped = _snv_alleles(variant, fasta_base)

            # Count each called allele. The site is biallelic, so a called allele is
            # 0 (REF) or 1 (ALT), and a missing one adds nothing.
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
            if swapped:
                swapped_at.add(pos_0)

    return snps, mismatched, swapped_at


def _query_outgroup_genotype(
    vcf: object,
    cds: CdsRegion,
    stats: GeneStats,
    ref_fetch: callable,
    include_filtered: bool,
) -> tuple[dict[int, tuple[str, ...]], set[int], set[int]]:
    """Query outgroup VCF for genotypes at CDS positions.

    Returns a dict mapping 0-based genomic position -> the base of each copy of the
    outgroup sample there, the positions whose REF is not the FASTA base, and the
    positions of calls read with REF and ALT swapped. The dict holds each call with every
    allele called and some ALT; each copy is an outgroup sequence, as in FASTA input. A
    position not in the dict is the reference. A missing or partly missing call is left
    out, treated like a missing record, and counted in ``stats.outgroup_skipped``.
    """
    outgroup_copies: dict[int, tuple[str, ...]] = {}
    mismatched: set[int] = set()
    swapped_at: set[int] = set()

    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"
        for variant in _query_region(vcf, region):
            pos_0 = variant.POS - 1

            if not _in_full_codon(cds, pos_0):
                continue

            fasta_base = ref_fetch(cds.chrom, pos_0)
            reason = _skip_reason(variant, fasta_base, include_filtered)
            if reason is _Skip.REF_MISMATCH:
                mismatched.add(pos_0)
            if reason is not None:
                if reason is not _Skip.INVARIANT:
                    stats.outgroup_skipped += 1
                continue
            alleles, ref, alt, swapped = _snv_alleles(variant, fasta_base)
            if alleles is None:
                stats.outgroup_skipped += 1
                continue

            # The outgroup is the first sample. Drop the padding past its ploidy.
            call = alleles[0][alleles[0] != -2]
            # A missing allele reads as REF, the same as no record; counted, as data is lost.
            if not (call.size and (call >= 0).all()):
                stats.outgroup_skipped += 1
                continue
            if (call == 1).any():
                outgroup_copies[pos_0] = tuple(alt if allele == 1 else ref for allele in call)
            if swapped:
                swapped_at.add(pos_0)

    return outgroup_copies, mismatched, swapped_at


def _divergence_outcome(
    in_codon: str, out_codons: set[str], clean_out: set[str], code: GeneticCode
) -> tuple[int, int] | _CodonSkip | None:
    """What one outgroup phase gives a codon: its (dn, ds), the counted reason it gives
    none, or None where the outgroup varies or matches the ingroup."""
    # An outgroup with more than one clean codon varies, so it has no fixed difference.
    if len(clean_out) > 1:
        return None
    # With no clean outgroup codon, any of them gives the skip reason below.
    out_codon = min(clean_out or out_codons)
    if in_codon == out_codon:
        return None
    # Checked before get_path, which has no path for such a codon and would otherwise
    # report the pair as stop-blocked.
    if not {in_codon, out_codon} <= CLEAN_CODONS:
        return _CodonSkip.AMBIGUOUS
    if {in_codon, out_codon} & code.stop_codons:
        return _CodonSkip.STOP
    # get_path handles codons that differ at more than one position.
    path = code.get_path(in_codon, out_codon)
    if not path:
        return _CodonSkip.STOP_BLOCKED
    kinds = [change_type for change_type, _position in path]
    return (kinds.count("R"), kinds.count("S"))


def _codon_with_bases(
    cds: CdsRegion, codon_idx: int, ref_fetch: callable, bases: dict[int, str]
) -> str:
    """The codon at ``codon_idx`` with ``bases`` in place of the reference at their positions."""
    return cds.extract_codon(codon_idx, lambda chrom, pos: bases.get(pos) or ref_fetch(chrom, pos))


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
    needs_derived_frequency: bool = False,
) -> tuple[PolymorphismData, GeneStats]:
    """Extract polymorphism and divergence data for a single gene from VCF.

    Args:
        vcf_path: Path to ingroup VCF file.
        outgroup_vcf_path: Path to outgroup VCF file, or None.
        cds: CDS region defining the gene.
        ref_fasta_path: Path to indexed reference FASTA.
        genetic_code: Genetic code for classification (standard if None).
        min_frequency: Minimum derived allele frequency threshold.
        no_singletons: If True, exclude each SNP whose derived allele is seen once.
        include_filtered: If True, also use records whose FILTER is not PASS or '.'.
        ingroup_vcf: Pre-opened cyvcf2.VCF handle for ingroup (optional).
        outgroup_vcf: Pre-opened cyvcf2.VCF handle for outgroup (optional).
        ref_fasta: Pre-opened pysam.FastaFile handle (optional).
        needs_derived_frequency: True when the analysis reads each SNP's derived
            frequency (the asymptotic, imputed and alpha-TG tests). A SNP the outgroup
            cannot polarize is then dropped and counted in ``unpolarized_sites``, as it
            is with ``min_frequency`` or ``no_singletons``; otherwise it is kept.

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
    ingroup_snps, mismatched, swapped_at = _query_ingroup_snps(
        vcf, cds, stats, ref_fetch, include_filtered
    )
    if ingroup_vcf is None:
        vcf.close()

    outgroup_copies: dict[int, tuple[str, ...]] = {}
    outgroup_swapped_at: set[int] = set()
    has_outgroup = outgroup_vcf is not None or outgroup_vcf_path is not None
    if has_outgroup:
        vcf = outgroup_vcf or _open_vcf(outgroup_vcf_path)
        outgroup_copies, outgroup_mismatched, outgroup_swapped_at = _query_outgroup_genotype(
            vcf, cds, stats, ref_fetch, include_filtered
        )
        mismatched |= outgroup_mismatched
        if outgroup_vcf is None:
            vcf.close()

    # The FASTA base at a REF mismatch cannot be trusted, so no record in a codon that
    # holds one is used, in either group.
    mismatch_codons = {cds.genomic_pos_to_codon_index(pos) for pos in mismatched}

    def usable(pos: int) -> bool:
        return cds.genomic_pos_to_codon_index(pos) not in mismatch_codons

    if mismatch_codons:
        kept_snps = [s for s in ingroup_snps if usable(s.pos)]
        stats.skipped_ref_mismatch += len(ingroup_snps) - len(kept_snps)
        ingroup_snps = kept_snps
        kept_copies = {pos: copies for pos, copies in outgroup_copies.items() if usable(pos)}
        stats.outgroup_skipped += len(outgroup_copies) - len(kept_copies)
        outgroup_copies = kept_copies
    # Counted after the codon rule, so a record is never both swapped and skipped.
    stats.swapped_ref_alt = sum(map(usable, swapped_at)) + sum(map(usable, outgroup_swapped_at))

    # === Polymorphism extraction ===
    polymorphisms: list[tuple[float, str]] = []
    unpolarized = 0
    needs_frequency = needs_derived_frequency or min_frequency > 0 or no_singletons
    # Every called sample carries a fixed ALT, so the background codon of the ingroup is the
    # reference with the fixed ALT bases, and each SNP is read against it, as in FASTA input.
    fixed_alt = {s.pos: s.alt for s in ingroup_snps if s.is_fixed_alt}

    def background_fetch(chrom: str, pos: int) -> str:
        return fixed_alt.get(pos) or ref_fetch(chrom, pos)

    def clean(codons: set[str]) -> set[str]:
        """The codons that are not missing data: no ambiguous base and not a stop."""
        return {c for c in codons if c in CLEAN_CODONS} - code.stop_codons

    # The outgroup records of each codon. Every record is in a full codon of the CDS.
    outgroup_positions: dict[int, list[int]] = defaultdict(list)
    for pos in sorted(outgroup_copies):
        outgroup_positions[cds.genomic_pos_to_codon_index(pos)].append(pos)

    def outgroup_phases(codon_idx: int) -> list[tuple[set[str], set[str]]]:
        """The outgroup's codons in each phase, all of them and the clean ones.

        Each copy of the outgroup sample gives one codon, as each outgroup sequence does in
        FASTA input. A call gives the base of each copy, but which copies of two calls are
        the same copy (the phase) is not known, so each way to pair them is a phase.
        """
        positions = outgroup_positions.get(codon_idx, [])
        ploidy = max((len(outgroup_copies[pos]) for pos in positions), default=1)
        built: dict[tuple[tuple[int, str], ...], str] = {}
        phases = []
        seen = set()
        # The copies have no names, so the first call's order is fixed, and only the others
        # are reordered against it.
        orders_per_position = [
            {outgroup_copies[pos]} if n == 0 else set(permutations(outgroup_copies[pos]))
            for n, pos in enumerate(positions)
        ]
        for orders in product(*orders_per_position):
            copies = []
            for copy in range(ploidy):
                bases = tuple(
                    (pos, order[copy % len(order)]) for pos, order in zip(positions, orders)
                )
                if bases not in built:
                    built[bases] = _codon_with_bases(cds, codon_idx, ref_fetch, dict(bases))
                copies.append(built[bases])
            # Phases that differ only in the order of the copies are the same phase.
            if (key := tuple(sorted(copies))) not in seen:
                seen.add(key)
                phases.append((set(copies), clean(set(copies))))
        return phases

    # Positions of SNPs dropped for a stop codon, in the ingroup or in the outgroup. Missing
    # data is not a polymorphism, so such a position does not void its codon for divergence.
    stop_allele_positions: set[int] = set()
    # SNPs whose background codon is a stop and whose ALT codon is not, alone among the
    # SNPs in their codon. The stop carriers are missing data, so the ALT base is the
    # ingroup's base there, as in FASTA input. With a second SNP in the codon, the samples
    # of each allele are not known, so no codon can be built and the stop stays missing data.
    alt_of_ref_stop: dict[int, str] = {}
    variants_per_codon = Counter(
        cds.genomic_pos_to_codon_index(s.pos) for s in ingroup_snps if not s.is_fixed_alt
    )
    # Where neither the outgroup nor a fixed ALT changes a codon, the outgroup codon is the
    # background codon, so only these codons need their own outgroup codon check. Without an
    # outgroup VCF there is no outgroup codon to check.
    outgroup_differs = (
        set(outgroup_positions) | {cds.genomic_pos_to_codon_index(pos) for pos in fixed_alt}
        if has_outgroup
        else set()
    )

    for snp in ingroup_snps:
        # Skip sites fixed for ALT in ingroup (these are potential divergence, not polymorphism)
        if snp.is_fixed_alt:
            continue

        codon_idx = cds.genomic_pos_to_codon_index(snp.pos)
        if codon_idx is None:
            continue

        # The codon without and with the SNP, on the ingroup background.
        ref_codon, alt_codon = _reconstruct_codon_with_sub(
            cds, codon_idx, background_fetch, snp.pos, snp.alt
        )

        # An ambiguous reference base (N, an IUPAC code) leaves the change unknown.
        if not {ref_codon, alt_codon} <= CLEAN_CODONS:
            stats.skipped_ambiguous += 1
            continue

        # A stop codon is missing data, so a change to or from one is not a polymorphism.
        if {ref_codon, alt_codon} & code.stop_codons:
            stats.skipped_stop += 1
            stop_allele_positions.add(snp.pos)
            if alt_codon not in code.stop_codons and variants_per_codon[codon_idx] == 1:
                alt_of_ref_stop[snp.pos] = snp.alt
            continue

        # An outgroup codon that is a stop or has an ambiguous base is missing data, as an
        # outgroup sequence with one is in FASTA input. With no clean outgroup codon in any
        # phase, the outgroup cannot polarize the polymorphism, and the FASTA path drops it.
        phases = [({ref_codon}, {ref_codon})]
        if codon_idx in outgroup_differs:
            phases = outgroup_phases(codon_idx)
            if not any(clean_codons for _, clean_codons in phases):
                all_out = set().union(*(codons for codons, _ in phases))
                if all_out & code.stop_codons:
                    # Divergence in this codon then meets the same stop codon and is skipped.
                    stats.skipped_stop += 1
                    stop_allele_positions.add(snp.pos)
                else:
                    stats.skipped_ambiguous += 1
                continue

        # Classify the change
        is_syn = code.is_synonymous_change(ref_codon, alt_codon)
        if is_syn is None:
            continue

        k = next(i for i in range(3) if ref_codon[i] != alt_codon[i])
        bases = {ref_codon[k]: snp.n_ref, alt_codon[k]: snp.n_alt}

        def outcome(clean_codons: set[str]) -> float | None | _SnpDrop:
            """For one outgroup phase: the derived frequency of a kept SNP, None for a kept
            SNP with no derived state, or why the SNP is dropped."""
            if not clean_codons:
                return _SnpDrop.MISSING
            # The ancestral base is the one the clean outgroup codons share, by the FASTA
            # rule. With none shared (a third allele) the SNP has no derived state: it
            # still counts in Pn/Ps, but not where a frequency is read.
            ancestral = ancestral_base(bases, clean_codons, k)
            if ancestral is None:
                return _SnpDrop.UNPOLARIZED if needs_frequency else None
            count = snp.n_ref if ancestral == alt_codon[k] else snp.n_alt
            state = (count / snp.total, count)
            if not passes_frequency_filter(state, min_frequency, no_singletons):
                return _SnpDrop.FILTERED
            return state[0]

        outcomes = [outcome(clean_codons) for _, clean_codons in phases]
        # The phase is not known, so an outcome that differs between phases is not used.
        # A kept SNP's frequency is part of its outcome only where the analysis reads it.
        keys = {
            result if isinstance(result, _SnpDrop) or needs_derived_frequency else None
            for result in outcomes
        }
        if len(keys) > 1:
            unpolarized += 1
            continue
        result = outcomes[0]
        if result is _SnpDrop.UNPOLARIZED:
            unpolarized += 1
        if isinstance(result, _SnpDrop):
            continue
        # Only the counts are read with no derived state, so any frequency will do.
        freq = snp.n_alt / snp.total if result is None else result
        polymorphisms.append((freq, "S" if is_syn else "N"))

    # === Divergence extraction ===
    dn = 0
    ds = 0

    if has_outgroup:
        # The ingroup codon carries ALT where the ingroup is fixed for it. A codon needs
        # one clean codon per group to be a fixed difference, so any polymorphic position,
        # including a third-allele outgroup base there, voids the whole codon, not just
        # that one base. A position with no outgroup record is the reference.
        polymorphic = {s.pos for s in ingroup_snps if not s.is_fixed_alt} - stop_allele_positions
        poly_codons = {cds.genomic_pos_to_codon_index(pos) for pos in polymorphic}

        def in_fetch(chrom: str, pos: int) -> str:
            return alt_of_ref_stop.get(pos) or background_fetch(chrom, pos)

        for codon_idx in sorted(outgroup_differs - poly_codons):
            in_codon = cds.extract_codon(codon_idx, in_fetch)
            outcomes = {
                _divergence_outcome(in_codon, codons, clean_codons, code)
                for codons, clean_codons in outgroup_phases(codon_idx)
            }
            # The phase is not known, so a fixed difference that differs between phases is
            # not used. A skip reason is counted only where every phase gives it.
            counts = {result if isinstance(result, tuple) else (0, 0) for result in outcomes}
            if len(counts) > 1:
                stats.skipped_unphased += 1
                continue
            ((codon_dn, codon_ds),) = counts
            dn += codon_dn
            ds += codon_ds
            if len(outcomes) == 1:
                match next(iter(outcomes)):
                    case _CodonSkip.AMBIGUOUS:
                        stats.skipped_ambiguous += 1
                    case _CodonSkip.STOP:
                        stats.skipped_stop += 1
                    case _CodonSkip.STOP_BLOCKED:
                        stats.stop_blocked_pairs += 1

    if owns_ref_fasta:
        ref_fasta.close()

    return PolymorphismData(
        polymorphisms=polymorphisms,
        dn=dn,
        ds=ds,
        gene_id=cds.gene_id,
        unpolarized_sites=unpolarized,
    ), stats
