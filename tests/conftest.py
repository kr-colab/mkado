"""Shared fixtures for the test suite.

The VCF fixtures build a small reference genome, indexed VCFs, and GFF3
annotations in temporary directories. The synthetic genome is designed so
that every gene has a known codon table and every variant record exercises
one specific parser branch, which lets tests assert exact counts.

Indexing uses pysam, so the bgzip and tabix binaries are not required.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pysam
import pytest

from mkado.core.cds import CdsRegion
from mkado.io.vcf import GeneStats

EXAMPLE_VCF_DIR = Path(__file__).parent.parent / "examples" / "example_vcf"

# chr1 holds five genes that between them cover plus and minus strands, a
# codon that spans an exon junction, a nonzero phase, a reference N, and a
# terminal stop codon. Spacers of T keep the genes apart.
CHR1_SEQ = (
    "ATGGCCAAATTCGGACTGAGCTACTAA"  # [0,27)  g_plus   +  ATG GCC AAA TTC GGA CTG AGC TAC TAA
    "TTT"
    "GAACTTGTCCAT"  # [30,42) g_minus  -  coding ATG GAC AAG TTC
    "TTT"
    "ATGAAAC"  # [45,52) g_split exon 1
    "GTAG"  # intron
    "CTGGG"  # [56,61) g_split exon 2; coding ATG AAA CCT GGG
    "TTTT"
    "AATGGTTCAC"  # [65,75) g_phase1 +, phase 1; coding ATG GTT CAC
    "TTT"
    "ATGCNTGGC"  # [78,87) g_ncodon +; coding ATG CNT GGC
    "TTT"
)
# chr2 holds twelve identical tiny genes so the CLI's parallel path (ten or more
# tasks) can run on synthetic data.
CHR2_SEQ = "ATGAAATTTCCC" * 12
# chrZ has VCF records but no FASTA entry. A gene on it makes the reference
# lookup fail inside a worker, which is the error path under test.
CONTIGS = {"chr1": len(CHR1_SEQ), "chr2": len(CHR2_SEQ), "chrZ": 50}

INGROUP_SAMPLES = ["s1", "s2", "s3", "s4"]
OUTGROUP_SAMPLES = ["outgroup"]

_VCF_PREAMBLE = [
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    *(f"##contig=<ID={name},length={length}>" for name, length in CONTIGS.items()),
]


def _write_vcf(
    path_stem: Path,
    records: list[str],
    *,
    samples: list[str] = INGROUP_SAMPLES,
    extra_header: list[str] = (),
) -> Path:
    """Write a VCF, then bgzip and tabix-index it.

    With no samples the FORMAT column is omitted, which gives a sites-only VCF.
    ``extra_header`` lines go after the standard preamble. The filename is
    built from ``path_stem.name`` so stems that contain a dot are not mangled.
    The uncompressed file stays beside the compressed one.
    """
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        columns += ["FORMAT", *samples]
    lines = [*_VCF_PREAMBLE, *extra_header, "\t".join(columns), *records]

    vcf_path = path_stem.parent / (path_stem.name + ".vcf")
    vcf_path.write_text("\n".join(lines) + "\n")
    gz_path = path_stem.parent / (path_stem.name + ".vcf.gz")
    pysam.tabix_compress(str(vcf_path), str(gz_path), force=True)
    pysam.tabix_index(str(gz_path), preset="vcf", force=True)
    return gz_path


@pytest.fixture
def write_vcf():
    """Expose the VCF writer to tests that need a one-off variant layout."""
    return _write_vcf


def _write_alignment(path: Path, records: dict[str, str]) -> Path:
    """Write a FASTA alignment from a name-to-sequence mapping."""
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
    return path


@pytest.fixture
def write_alignment():
    return _write_alignment


# ---------------------------------------------------------------------------
# Synthetic genome
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GeneSpec:
    """Geometry of one synthetic gene. Exons are 0-based half-open, genomic order."""

    gene_id: str
    chrom: str
    exons: tuple[tuple[int, int], ...]
    strand: str = "+"
    phase: int = 0

    @property
    def transcript_id(self) -> str:
        return f"tx_{self.gene_id}"

    def cds(self) -> CdsRegion:
        return CdsRegion(
            gene_id=self.gene_id,
            transcript_id=self.transcript_id,
            chrom=self.chrom,
            exons=list(self.exons),
            strand=self.strand,
            phase=self.phase,
        )


def _gff_line(spec: GeneSpec, kind: str, start: int, end: int, phase: object, attrs: str) -> str:
    return "\t".join(
        [spec.chrom, ".", kind, str(start), str(end), ".", spec.strand, str(phase), attrs]
    )


def _write_gff3(path: Path, specs: tuple[GeneSpec, ...]) -> Path:
    """Write a GFF3 with one gene, mRNA, and CDS block per spec.

    CDS lines are emitted in coding order and carry the GFF phase of each exon,
    which is the number of bases to skip before the first complete codon.
    """
    lines = ["##gff-version 3"]
    for spec in specs:
        gene_start = min(start for start, _ in spec.exons) + 1
        gene_end = max(end for _, end in spec.exons)
        gene_attrs = f"ID={spec.gene_id};Name={spec.gene_id}"
        mrna_attrs = f"ID={spec.transcript_id};Parent={spec.gene_id}"
        lines.append(_gff_line(spec, "gene", gene_start, gene_end, ".", gene_attrs))
        lines.append(_gff_line(spec, "mRNA", gene_start, gene_end, ".", mrna_attrs))
        phase = spec.phase
        for i, (start, end) in enumerate(sorted(spec.exons, reverse=spec.strand == "-")):
            attrs = f"ID=cds_{spec.gene_id}_{i};Parent={spec.transcript_id}"
            lines.append(_gff_line(spec, "CDS", start + 1, end, phase, attrs))
            phase = (3 - ((end - start) - phase) % 3) % 3
    path.write_text("\n".join(lines) + "\n")
    return path


CHR1_GENES = (
    GeneSpec("g_plus", "chr1", ((0, 27),)),
    GeneSpec("g_minus", "chr1", ((30, 42),), strand="-"),
    GeneSpec("g_split", "chr1", ((45, 52), (56, 61))),
    GeneSpec("g_phase1", "chr1", ((65, 75),), phase=1),
    GeneSpec("g_ncodon", "chr1", ((78, 87),)),
)
TINY_GENES = tuple(GeneSpec(f"tiny{i:02d}", "chr2", ((12 * i, 12 * i + 9),)) for i in range(12))
ALL_GENES = CHR1_GENES + TINY_GENES
BAD_CHROM_GENE = GeneSpec("g_nochrom", "chrZ", ((0, 9),))


def _rec(chrom: str, pos: int, ref: str, alt: str, gts: list[str], info: str = ".") -> str:
    return "\t".join([chrom, str(pos), ".", ref, alt, "30", "PASS", info, "GT", *gts])


def _chr2_ingroup_records() -> list[str]:
    """Every tiny gene gets a synonymous singleton; even ones also get a nonsynonymous site."""
    records = []
    for i in range(12):
        base = 12 * i
        if i % 2 == 0:
            records.append(
                _rec("chr2", base + 8, "T", "A", ["0/1", "0/1", "0/0", "0/0"])
            )  # TTT>TAT
        records.append(_rec("chr2", base + 9, "T", "C", ["0/1", "0/0", "0/0", "0/0"]))  # TTT>TTC
    return records


# Each record targets one parser branch. Positions are 1-based. Records stay
# sorted by contig block and position because tabix requires it.
INGROUP_RECORDS = [
    _rec("chr1", 1, "AT", "A", ["0/1", "0/0", "0/0", "0/0"]),  # indel, skipped
    _rec("chr1", 3, "G", "A,T", ["0/1", "0/2", "0/0", "0/0"]),  # multi-allelic, skipped
    _rec("chr1", 6, "C", "T", ["0/0", "0/0", "0/0", "0/1"]),  # GCC>GCT  S  0.125
    _rec("chr1", 8, "A", "G", ["0/0", "0/0", "0/1", "0/1"]),  # AAA>AGA  N  0.25
    _rec("chr1", 10, "T", "<DEL>", ["0/1", "0/0", "0/0", "0/0"], "SVTYPE=DEL"),  # symbolic
    _rec("chr1", 11, "T", "C", ["1/1"] * 4),  # fixed for ALT: TTC>TCC against the outgroup
    _rec("chr1", 13, "G", "T", ["./."] * 4),  # all genotypes missing
    _rec("chr1", 17, "T", "C", ["0/0", "0/0", "0/1", "0/1"]),  # CTG>CCG  N, outgroup has C
    _rec("chr1", 22, "T", "C", ["0/0"] * 4),  # ALT frequency zero
    _rec("chr1", 27, "A", "G", ["0/0", "0/0", "0/0", "0/1"]),  # inside reference stop codon
    _rec("chr1", 34, "C", "T", ["0/1", "0/0", "0/0", "0/0"]),  # g_minus AAG>AAA  S  0.125
    _rec("chr1", 38, "T", "G", ["0/0", "0/1", "0/1", "0/1"]),  # g_minus GAC>GCC  N  0.375
    _rec("chr1", 54, "T", "C", ["0/1"] * 4),  # g_split intron, never queried
    _rec("chr1", 57, "C", "A", ["0/0", "0/1", "0/0", "0/0"]),  # g_split CCT>CAT  N  0.125
    _rec("chr1", 66, "A", "G", ["0/1", "0/0", "0/0", "0/0"]),  # g_phase1 trimmed base
    _rec("chr1", 72, "T", "C", ["0/1", "0/1", "0/0", "0/0"]),  # g_phase1 GTT>GTC  S  0.25
    _rec("chr1", 84, "T", "A", ["0/1", "0/0", "0/0", "0/0"]),  # g_ncodon CNT>CNA, unclassifiable
    _rec("chr1", 87, "C", "T", ["0/1", "0/1", "0/1", "0/0"]),  # g_ncodon GGC>GGT  S  0.375
    *_chr2_ingroup_records(),
    _rec("chrZ", 5, "A", "G", ["0/1", "0/0", "0/0", "0/0"]),  # contig absent from the FASTA
]

OUTGROUP_RECORDS = [
    _rec("chr1", 5, "CC", "C", ["1/1"]),  # indel, ignored
    _rec("chr1", 7, "A", "G,T", ["1/1"]),  # multi-allelic, ignored
    _rec("chr1", 9, "A", "<DEL>", ["1/1"], "SVTYPE=DEL"),  # symbolic, ignored
    _rec("chr1", 12, "C", "T", ["1/1"]),  # TTC>TTT  ds
    _rec("chr1", 14, "G", "A", ["0/1"]),  # GGA>GAA, heterozygous outgroup, discarded
    _rec("chr1", 17, "T", "C", ["1/1"]),  # shared with an ingroup polymorphism, discarded
    _rec("chr1", 18, "G", "A", ["1/1"]),  # CTG>CTA  ds
    _rec("chr1", 19, "A", "G", ["0/0"]),  # homozygous reference, ignored
    _rec("chr1", 21, "T", "C", ["1/1"]),  # ALT equals the FASTA base, codons identical
    _rec("chr1", 24, "C", "A", ["1/1"]),  # TAC>TAA creates a stop, ignored
    _rec("chr1", 26, "A", "G", ["1/1"]),  # TAA>TGA inside the reference stop, ignored
    _rec("chr1", 31, "G", "A", ["1/1"]),  # g_minus TTC>TTT  ds
    _rec("chr1", 40, "C", "T", ["1/1"]),  # g_minus ATG>ATA  dn (table 2: ds)
    _rec(
        "chr1", 52, "C", "T", ["1/1"]
    ),  # g_split CCT>TCT, shares codon with ingroup poly, discarded
    _rec("chr1", 61, "G", "A", ["1/1"]),  # g_split GGG>GGA  ds
    _rec("chr1", 66, "A", "G", ["1/1"]),  # g_phase1 trimmed base, ignored
    _rec("chr1", 74, "A", "G", ["1/1"]),  # g_phase1 CAC>CGC  dn
    _rec("chr1", 82, "C", "A", ["1/1"]),  # g_ncodon CNT>ANT, no path through X
    *(_rec("chr2", 12 * i + 5, "A", "G", ["1/1"]) for i in range(12)),  # AAA>AGA  dn
]


@dataclass(frozen=True)
class Expected:
    """What extract_gene_data returns for one gene of the synthetic genome."""

    polymorphisms: tuple[tuple[float, str], ...]
    dn: int
    ds: int
    skipped: tuple[int, int, int] = (0, 0, 0)  # indels, multi-allelic, missing
    warning: str | None = None

    @property
    def stats(self) -> GeneStats:
        return GeneStats(*self.skipped)

    @property
    def counts(self) -> tuple[int, int, int, int]:
        """(dn, ds, pn, ps)"""
        pn = sum(1 for _, kind in self.polymorphisms if kind == "N")
        return (self.dn, self.ds, pn, len(self.polymorphisms) - pn)


EXPECTED = {
    "g_plus": Expected(
        ((0.125, "S"), (0.25, "N"), (0.75, "N")),
        dn=1,
        ds=1,
        skipped=(2, 1, 1),
        warning="g_plus: skipped 2 indels, 1 multi-allelic",
    ),
    "g_minus": Expected(((0.125, "S"), (0.375, "N")), dn=1, ds=1),
    "g_split": Expected(((0.125, "N"),), dn=0, ds=1),
    "g_phase1": Expected(((0.25, "S"),), dn=1, ds=0),
    "g_ncodon": Expected(((0.375, "S"),), dn=0, ds=0),
}


@dataclass(frozen=True)
class VcfDataset:
    """One ingroup, outgroup, reference triple plus what the tests know about it.

    Tests override a path with ``dataclasses.replace``.
    """

    ref_fasta: Path
    ingroup_vcf: Path
    outgroup_vcf: Path
    gff: Path | None = None
    genes: dict[str, GeneSpec] = field(default_factory=dict)
    expected: dict[str, Expected] = field(default_factory=dict)

    def cds(self, gene_id: str) -> CdsRegion:
        return self.genes[gene_id].cds()


@pytest.fixture(scope="session")
def genome(tmp_path_factory: pytest.TempPathFactory) -> VcfDataset:
    """The synthetic genome, built once per session and treated as read-only."""
    root = tmp_path_factory.mktemp("genome")
    ref = root / "ref.fa"
    ref.write_text(f">chr1\n{CHR1_SEQ}\n>chr2\n{CHR2_SEQ}\n")
    pysam.faidx(str(ref))
    ingroup = _write_vcf(root / "ingroup", INGROUP_RECORDS)
    outgroup = _write_vcf(root / "outgroup", OUTGROUP_RECORDS, samples=OUTGROUP_SAMPLES)
    genes = {spec.gene_id: spec for spec in (*ALL_GENES, BAD_CHROM_GENE)}
    return VcfDataset(ref, ingroup, outgroup, genes=genes, expected=EXPECTED)


@pytest.fixture
def ingroup_vcf_plain(genome: VcfDataset) -> Path:
    """The uncompressed, unindexed ingroup VCF that the writer leaves beside the indexed one."""
    return genome.ingroup_vcf.with_suffix("")


@pytest.fixture
def outgroup_vcf_plain(genome: VcfDataset) -> Path:
    return genome.outgroup_vcf.with_suffix("")


@pytest.fixture
def outgroup_vcf_sites_only(tmp_path: Path) -> Path:
    """An outgroup VCF with no FORMAT or sample columns."""
    return _write_vcf(tmp_path / "outgroup_sites", ["chr1\t12\t.\tC\tT\t30\tPASS\t."], samples=[])


@pytest.fixture
def not_a_vcf(tmp_path: Path) -> Path:
    path = tmp_path / "not.vcf"
    path.write_text("hello\n")
    return path


@pytest.fixture
def bgzipped_ref(tmp_path: Path, genome: VcfDataset) -> Path:
    """The synthetic reference as a bgzipped, faidx-indexed FASTA."""
    gz = tmp_path / "ref.fa.gz"
    pysam.tabix_compress(str(genome.ref_fasta), str(gz), force=True)
    pysam.faidx(str(gz))
    return gz


@pytest.fixture
def gff_chr1(tmp_path: Path) -> Path:
    return _write_gff3(tmp_path / "chr1.gff3", CHR1_GENES)


@pytest.fixture
def gff_all(tmp_path: Path) -> Path:
    return _write_gff3(tmp_path / "all.gff3", ALL_GENES)


@pytest.fixture
def gff_tiny(tmp_path: Path) -> Path:
    return _write_gff3(tmp_path / "tiny.gff3", TINY_GENES)


@pytest.fixture
def gff_bad_chrom(tmp_path: Path) -> Path:
    return _write_gff3(tmp_path / "bad_chrom.gff3", (BAD_CHROM_GENE, CHR1_GENES[0]))


@pytest.fixture
def gff_all_with_bad_chrom(tmp_path: Path) -> Path:
    """Enough genes for the parallel path, with one that fails inside its chunk."""
    return _write_gff3(tmp_path / "all_bad_chrom.gff3", (*ALL_GENES, BAD_CHROM_GENE))


@pytest.fixture
def gff_bad_chrom_only(tmp_path: Path) -> Path:
    return _write_gff3(tmp_path / "bad_chrom_only.gff3", (BAD_CHROM_GENE,))


@pytest.fixture
def gff_invalid_len(tmp_path: Path) -> Path:
    """A GFF3 whose only gene has a CDS length that is not a multiple of three."""
    return _write_gff3(tmp_path / "invalid.gff3", (GeneSpec("g_bad", "chr1", ((0, 10),)),))


@pytest.fixture
def fitter_calls(monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    """Record the keyword arguments of each asymptotic fit made in this process.

    Covers both the aggregated fitter (used by the VCF path for every gene,
    and by the FASTA aggregated batch mode) and the FASTA per-gene fitter.
    The workers import the fitter lazily, so patching the module attribute
    reaches them.
    """
    import mkado.analysis.asymptotic as asymptotic

    calls: list[dict] = []

    def _recorder(real):
        def recording(*args, **kwargs):
            calls.append(kwargs)
            return real(*args, **kwargs)

        return recording

    for name in ("asymptotic_mk_test_aggregated", "asymptotic_mk_test"):
        monkeypatch.setattr(asymptotic, name, _recorder(getattr(asymptotic, name)))
    return calls


@pytest.fixture(scope="session")
def example_dataset() -> VcfDataset:
    """The tracked example data: 42 genes, bgzipped reference, indexed VCFs."""
    if not EXAMPLE_VCF_DIR.exists():
        pytest.skip("example VCF data not available")
    return VcfDataset(
        ref_fasta=EXAMPLE_VCF_DIR / "reference.fa.gz",
        ingroup_vcf=EXAMPLE_VCF_DIR / "ingroup.vcf.gz",
        outgroup_vcf=EXAMPLE_VCF_DIR / "outgroup.vcf.gz",
        gff=EXAMPLE_VCF_DIR / "annotation.gff3",
    )
