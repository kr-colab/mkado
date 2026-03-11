"""VCF polymorphism and divergence extraction for MK tests."""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from mkado.analysis.asymptotic import PolymorphismData
from mkado.core.cds import CdsRegion, _COMPLEMENT
from mkado.core.codons import DEFAULT_CODE, GeneticCode

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def _check_vcf_deps() -> None:
    """Check that cyvcf2 and pysam are available."""
    try:
        import cyvcf2  # noqa: F401
    except ImportError:
        raise ImportError("cyvcf2 is required for VCF input. Install with: pip install mkado[vcf]")
    try:
        import pysam  # noqa: F401
    except ImportError:
        raise ImportError("pysam is required for VCF input. Install with: pip install mkado[vcf]")


@dataclass
class _SnpInfo:
    """A biallelic SNP in a CDS region."""

    pos: int  # 0-based genomic position
    ref: str
    alt: str
    alt_freq: float  # ALT allele frequency in ingroup
    n_samples: int  # number of non-missing samples
    is_fixed_alt: bool  # all ingroup samples are ALT


@dataclass
class GeneStats:
    """Bookkeeping stats for a single gene's VCF extraction."""

    skipped_indels: int = 0
    skipped_multiallelic: int = 0
    skipped_missing: int = 0


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
    vcf_path: str | Path,
    cds: CdsRegion,
    stats: GeneStats,
) -> list[_SnpInfo]:
    """Query ingroup VCF for biallelic SNPs overlapping a CDS region.

    Returns SNPs in genomic coordinates with allele frequencies.
    """
    import cyvcf2

    vcf = cyvcf2.VCF(str(vcf_path))
    snps: list[_SnpInfo] = []

    # Build query regions from exons
    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"  # cyvcf2 uses 1-based
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                variants = list(vcf(region))
        except Exception:
            continue

        for variant in variants:
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

            alt_freq = n_alt / total
            n_samples = sum(1 for gt in gt_types if gt != 2)
            is_fixed_alt = n_ref == 0

            snps.append(
                _SnpInfo(
                    pos=pos_0,
                    ref=ref,
                    alt=alt,
                    alt_freq=alt_freq,
                    n_samples=n_samples,
                    is_fixed_alt=is_fixed_alt,
                )
            )

    vcf.close()
    return snps


def _query_outgroup_genotype(
    vcf_path: str | Path,
    cds: CdsRegion,
) -> dict[int, str]:
    """Query outgroup VCF for genotypes at CDS positions.

    Returns dict mapping 0-based genomic position -> outgroup allele.
    Only includes positions where outgroup differs from reference.
    """
    import cyvcf2

    vcf = cyvcf2.VCF(str(vcf_path))
    outgroup_alleles: dict[int, str] = {}

    for start, end in cds.exons:
        region = f"{cds.chrom}:{start + 1}-{end}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                variants = list(vcf(region))
        except Exception:
            continue

        for variant in variants:
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
            elif gt == 1:  # HET - treat as carrying ALT
                outgroup_alleles[pos_0] = alt

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
) -> tuple[PolymorphismData, GeneStats]:
    """Extract polymorphism and divergence data for a single gene from VCF.

    Args:
        vcf_path: Path to ingroup VCF file.
        outgroup_vcf_path: Path to outgroup VCF file, or None.
        cds: CDS region defining the gene.
        ref_fasta_path: Path to indexed reference FASTA.
        genetic_code: Genetic code for classification (standard if None).
        min_frequency: Minimum derived allele frequency threshold.
        no_singletons: If True, exclude singletons.

    Returns:
        Tuple of (PolymorphismData, GeneStats).
    """
    import pysam

    code = genetic_code or DEFAULT_CODE
    stats = GeneStats()

    ref_fasta = pysam.FastaFile(str(ref_fasta_path))
    ref_fetch = _ref_base_fetcher(ref_fasta)

    # Get ingroup SNPs
    ingroup_snps = _query_ingroup_snps(vcf_path, cds, stats)

    # Calculate singleton threshold
    if no_singletons and ingroup_snps:
        max_n = max(s.n_samples for s in ingroup_snps)
        singleton_freq = 1.0 / (2 * max_n)  # diploid
        if singleton_freq > min_frequency:
            min_frequency = singleton_freq

    # Get outgroup genotypes
    outgroup_alleles: dict[int, str] = {}
    if outgroup_vcf_path is not None:
        outgroup_alleles = _query_outgroup_genotype(outgroup_vcf_path, cds)

    # Index ingroup SNPs by position for quick lookup
    ingroup_by_pos: dict[int, _SnpInfo] = {s.pos: s for s in ingroup_snps}

    # Track which codon positions have ingroup variants (for divergence check)
    ingroup_variant_positions: set[int] = {s.pos for s in ingroup_snps}

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

        # Determine derived allele frequency
        # If outgroup carries ALT, then REF is derived
        derived_freq = snp.alt_freq
        if snp.pos in outgroup_alleles and outgroup_alleles[snp.pos] == snp.alt:
            derived_freq = 1.0 - snp.alt_freq

        # Apply frequency filter
        if derived_freq < min_frequency:
            continue
        if derived_freq <= 0.0 or derived_freq >= 1.0:
            continue

        poly_type = "S" if is_syn else "N"
        polymorphisms.append((derived_freq, poly_type))

    # === Divergence extraction ===
    dn = 0
    ds = 0

    if outgroup_vcf_path is not None:
        # Check each codon for fixed differences between reference and outgroup
        for codon_idx in range(cds.num_codons()):
            p1, p2, p3 = cds.codon_positions(codon_idx)

            # Check if outgroup differs at any position in this codon
            has_outgroup_diff = False
            for pos in (p1, p2, p3):
                if pos in outgroup_alleles:
                    # Only count as divergence if ingroup is monomorphic for REF at this position
                    if pos in ingroup_variant_positions:
                        ingroup_snp = ingroup_by_pos.get(pos)
                        if ingroup_snp is not None and not ingroup_snp.is_fixed_alt:
                            # Ingroup is polymorphic — not a fixed divergence
                            continue
                    has_outgroup_diff = True

            if not has_outgroup_diff:
                continue

            # Reconstruct reference and outgroup codons
            ref_bases = [ref_fetch(cds.chrom, p) for p in (p1, p2, p3)]
            out_bases = list(ref_bases)

            for i, pos in enumerate((p1, p2, p3)):
                if pos in outgroup_alleles:
                    # Check ingroup is monomorphic for REF at this position
                    if pos in ingroup_variant_positions:
                        ingroup_snp = ingroup_by_pos.get(pos)
                        if ingroup_snp is not None and not ingroup_snp.is_fixed_alt:
                            out_bases[i] = ref_bases[i]  # Don't count this position
                            continue
                    out_bases[i] = outgroup_alleles[pos]

            ref_codon = "".join(ref_bases)
            out_codon = "".join(out_bases)

            if cds.strand == "-":
                ref_codon = ref_codon.translate(_COMPLEMENT)
                out_codon = out_codon.translate(_COMPLEMENT)

            ref_codon = ref_codon.upper()
            out_codon = out_codon.upper()

            if ref_codon == out_codon:
                continue

            # Skip stop codons
            if code.translate(ref_codon) == "*" or code.translate(out_codon) == "*":
                continue

            # Classify using get_path for multi-step changes
            path = code.get_path(ref_codon, out_codon)
            for change_type, _position in path:
                if change_type == "R":
                    dn += 1
                elif change_type == "S":
                    ds += 1

    ref_fasta.close()

    return PolymorphismData(
        polymorphisms=polymorphisms,
        dn=dn,
        ds=ds,
        gene_id=cds.gene_id,
    ), stats
