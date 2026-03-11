# MKado 御門

[![Documentation Status](https://readthedocs.org/projects/mkado/badge/?version=latest)](https://mkado.readthedocs.io/en/latest/?badge=latest)

A modern Python implementation of the McDonald-Kreitman test toolkit for detecting selection in molecular evolution.

**[Documentation](https://mkado.readthedocs.io/)** | **[PyPI](https://pypi.org/project/mkado/)**

## Features

- **Standard MK test**: Classic 2x2 contingency table with Fisher's exact test
- **Polarized MK test**: Uses a third outgroup to assign mutations to lineages
- **Asymptotic MK test**: Frequency-bin α estimates with exponential extrapolation (Messer & Petrov 2013)
- **Tarone-Greenland α_TG**: Weighted multi-gene estimator that corrects for sample size heterogeneity (Stoletzki & Eyre-Walker 2011)
- **Alternate genetic codes**: 24 NCBI genetic code tables (mitochondrial, plastid, etc.) selectable by name
- **VCF input**: Go directly from VCF + reference + GFF3 annotation to MK test results (no FASTA alignment needed)
- **Batch processing**: Process multiple genes with parallel execution and Benjamini-Hochberg correction for multiple testing
- **Volcano plots**: Visualize batch results with publication-ready volcano plots
- **Multiple output formats**: Pretty-print, TSV, and JSON

## Installation

```bash
# Install with uv (recommended)
uv pip install mkado

# Or install with pip
pip install mkado

# With VCF support (adds cyvcf2 and pysam)
pip install mkado[vcf]
```

### Development Installation

```bash
# Clone the repository
git clone https://github.com/andrewkern/mkado.git
cd mkado

# Install with uv
uv sync
```

## Quick Start

```bash
# Standard MK test (combined alignment file)
mkado test alignment.fa -i "dmel" -o "dsim"

# Asymptotic MK test
mkado test alignment.fa -i "dmel" -o "dsim" -a

# Polarized MK test
mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"

# Batch process a directory
mkado batch alignments/ -i "dmel" -o "dsim"

# Batch with asymptotic test and 8 parallel workers
mkado batch alignments/ -i "dmel" -o "dsim" -a -w 8

# Run MK test directly from VCF files
mkado vcf --vcf pop.vcf.gz --outgroup-vcf out.vcf.gz --ref genome.fa --gff genes.gff3

# Get file info
mkado info sequences.fa
```

## Usage Modes

mkado supports two modes for specifying ingroup/outgroup sequences:

### Combined File Mode (Recommended)

Use `-i` and `-o` to filter sequences by name pattern from a single alignment file:

```bash
mkado test alignment.fa -i "speciesA" -o "speciesB"
mkado batch alignments/ -i "speciesA" -o "speciesB"
```

### Separate Files Mode

Provide separate FASTA files for ingroup and outgroup:

```bash
mkado test ingroup.fa outgroup.fa
mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"
```

## Commands

### `mkado test`

Run MK test on a single alignment.

```bash
mkado test FASTA [OUTGROUP_FILE] [OPTIONS]
```

**Key Options:**
| Option | Short | Description |
|--------|-------|-------------|
| `--ingroup-match` | `-i` | Ingroup sequence name pattern (combined mode) |
| `--outgroup-match` | `-o` | Outgroup sequence name pattern (combined mode) |
| `--asymptotic` | `-a` | Use asymptotic MK test |
| `--polarize` | `-p` | Second outgroup file (separate files mode) |
| `--polarize-match` | | Second outgroup pattern (combined mode) |
| `--bins` | `-b` | Frequency bins for asymptotic test (default: 10) |
| `--plot-asymptotic` | | Generate alpha(x) plot for asymptotic test (PNG, PDF, or SVG) |
| `--code-table` | | Genetic code (e.g. `vertebrate-mito`, or NCBI ID) |
| `--format` | `-f` | Output format: pretty, tsv, json |
| `--reading-frame` | `-r` | Reading frame 1-3 (default: 1) |

**Examples:**

```bash
# Combined file mode
mkado test alignment.fa -i "dmel" -o "dsim"
mkado test alignment.fa -i "dmel" -o "dsim" -a -b 20
mkado test alignment.fa -i "dmel" -o "dsim" -a --plot-asymptotic alpha_fit.png
mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"

# Separate files mode
mkado test ingroup.fa outgroup.fa
mkado test ingroup.fa outgroup.fa -a
mkado test ingroup.fa outgroup.fa -p outgroup2.fa
```

### `mkado batch`

Run MK test on multiple alignment files.

```bash
mkado batch DIRECTORY [OPTIONS]
```

**Key Options:**
| Option | Short | Description |
|--------|-------|-------------|
| `--ingroup-match` | `-i` | Ingroup pattern (enables combined file mode) |
| `--outgroup-match` | `-o` | Outgroup pattern (required with -i) |
| `--asymptotic` | `-a` | Use asymptotic MK test |
| `--alpha-tg` | | Compute weighted α_TG (Stoletzki & Eyre-Walker 2011) |
| `--aggregate/--per-gene` | | Aggregate results or per-gene (asymptotic) |
| `--pattern` | | File glob pattern (default: auto-detect *.fa, *.fasta, *.fna) |
| `--workers` | `-w` | Parallel workers (0=auto, 1=sequential) |
| `--bins` | `-b` | Frequency bins for asymptotic test |
| `--code-table` | | Genetic code (e.g. `vertebrate-mito`, or NCBI ID) |
| `--format` | `-f` | Output format: pretty, tsv, json |
| `--volcano` | | Generate volcano plot (PNG, PDF, or SVG) |
| `--plot-asymptotic` | | Generate alpha(x) plot for aggregated asymptotic test |

**Examples:**

```bash
# Combined file mode (recommended)
mkado batch alignments/ -i "dmel" -o "dsim"
mkado batch alignments/ -i "dmel" -o "dsim" -a
mkado batch alignments/ -i "dmel" -o "dsim" -a --per-gene
mkado batch alignments/ -i "dmel" -o "dsim" --alpha-tg
mkado batch alignments/ -i "dmel" -o "dsim" -w 8

# Generate a volcano plot
mkado batch alignments/ -i "dmel" -o "dsim" --volcano results.png

# Generate asymptotic alpha(x) plot
mkado batch alignments/ -i "dmel" -o "dsim" -a --plot-asymptotic alpha_fit.png

# Separate files mode
mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"
```

### `mkado vcf`

Run MK test from VCF + reference genome + GFF3 annotation. No pre-aligned FASTA files needed.

```bash
mkado vcf --vcf VCF --outgroup-vcf VCF --ref FASTA --gff GFF3 [OPTIONS]
```

**Key Options:**
| Option | Short | Description |
|--------|-------|-------------|
| `--vcf` | | Ingroup VCF (multi-sample, bgzipped+tabix recommended) |
| `--outgroup-vcf` | | Single-sample outgroup VCF |
| `--ref` | | Reference FASTA (plain or bgzipped, faidx-indexed) |
| `--gff` | | GFF3 annotation (plain or gzipped) |
| `--gene` | | Analyze a single gene by ID |
| `--gene-list` | | File with gene IDs to analyze |
| `--asymptotic` | `-a` | Use asymptotic MK test |
| `--alpha-tg` | | Compute weighted α_TG |
| `--imputed` | | Use imputed MK test |
| `--aggregate/--per-gene` | | Aggregate or per-gene results |
| `--workers` | `-w` | Parallel workers (0=auto) |
| `--format` | `-f` | Output format: pretty, tsv, json |

**Examples:**

```bash
# Standard MK test across all genes
mkado vcf --vcf pop.vcf.gz --outgroup-vcf out.vcf.gz --ref genome.fa --gff genes.gff3

# Asymptotic MK test
mkado vcf --vcf pop.vcf.gz --outgroup-vcf out.vcf.gz --ref genome.fa --gff genes.gff3 -a

# Single gene analysis
mkado vcf --vcf pop.vcf.gz --outgroup-vcf out.vcf.gz --ref genome.fa --gff genes.gff3 --gene BRCA1

# With parallel workers and TSV output
mkado vcf --vcf pop.vcf.gz --outgroup-vcf out.vcf.gz --ref genome.fa --gff genes.gff3 -w 8 -f tsv
```

### `mkado codes`

List available genetic code tables.

```bash
mkado codes
```

### `mkado info`

Display information about a FASTA file.

```bash
mkado info FASTA [-r READING_FRAME]
```

## Example Output

```
$ mkado test alignment.fa -i "kreitman" -o "mauritiana"

Found 11 ingroup, 1 outgroup sequences
MK Test Results:
  Divergence:    Dn=6, Ds=8
  Polymorphism:  Pn=1, Ps=8
  Fisher's exact p-value: 0.176
  Neutrality Index (NI):  0.1667
  Alpha (α):              0.8333
  DoS:                    0.3175
```

## Python API

```python
from mkado import mk_test, asymptotic_mk_test, SequenceSet

# Run MK test
result = mk_test("ingroup.fa", "outgroup.fa")
print(f"Alpha: {result.alpha}")
print(f"P-value: {result.p_value}")

# Run asymptotic MK test
result = asymptotic_mk_test("ingroup.fa", "outgroup.fa")
print(f"Asymptotic Alpha: {result.alpha_asymptotic}")
print(f"95% CI: {result.ci_low} - {result.ci_high}")

# Combined file mode - filter by sequence name
all_seqs = SequenceSet.from_fasta("combined.fa")
ingroup = all_seqs.filter_by_name("dmel")
outgroup = all_seqs.filter_by_name("dsim")
result = mk_test(ingroup, outgroup)
```

## Interpretation

### Neutrality Index (NI)

- **NI = 1**: Neutral evolution
- **NI > 1**: Excess polymorphism (segregating weakly deleterious variants)
- **NI < 1**: Excess divergence (positive selection)

### Alpha (α)

- **α = 0**: No adaptive substitutions
- **α > 0**: Proportion of substitutions driven by positive selection
- **α < 0**: Excess polymorphism relative to divergence

### Direction of Selection (DoS)

From Stoletzki & Eyre-Walker (2011), DoS = Dn/(Dn+Ds) - Pn/(Pn+Ps):

- **DoS = 0**: Neutral evolution
- **DoS > 0**: Positive selection (excess adaptive substitutions)
- **DoS < 0**: Slightly deleterious polymorphisms

DoS is bounded [-1, +1] and symmetric around 0, making it easier to interpret than NI.

## Development

```bash
# Install dev dependencies
uv sync

# Run tests
uv run pytest

# Run linter
uv run ruff check src/

# Run formatter
uv run ruff format src/
```

## Examples

Example data and tutorials are available in the `examples/` directory:

```bash
# Run batch MK test on example data
mkado batch examples/anopheles_batch/ -i gamb -o afun

# Run asymptotic MK test
mkado batch examples/anopheles_batch/ -i gamb -o afun -a
```

See the [documentation](https://mkado.readthedocs.io/) for detailed tutorials and API reference.

## Citation

If you use MKado in your research, please cite:

> Rivera-Colón, A. G., Rehmann, C. T., & Kern, A. D. (2026). MKado: a toolkit for McDonald-Kreitman tests of natural selection. *bioRxiv*. https://doi.org/10.64898/2026.03.02.709122

## References

- Haller, B. C., & Messer, P. W. (2017). asymptoticMK: A web-based tool for the asymptotic McDonald–Kreitman test. G3: Genes, Genomes, Genetics, 7(5), 1569-1575. https://doi.org/10.1534/g3.117.039693
- McDonald, J. H., & Kreitman, M. (1991). Adaptive protein evolution at the Adh locus in Drosophila. Nature, 351(6328), 652-654. https://doi.org/10.1038/351652a0
- Messer, P. W., & Petrov, D. A. (2013). Frequent adaptation and the McDonald–Kreitman test. PNAS, 110(21), 8615-8620. https://doi.org/10.1073/pnas.1220835110
- Rand, D. M., & Kann, L. M. (1996). Excess amino acid polymorphism in mitochondrial DNA: contrasts among genes from Drosophila, mice, and humans. Molecular Biology and Evolution, 13(6), 735-748. https://doi.org/10.1093/oxfordjournals.molbev.a025634
- Smith, N. G. C., & Eyre-Walker, A. (2002). Adaptive protein evolution in Drosophila. Nature, 415(6875), 1022-1024. https://doi.org/10.1038/4151022a
- Stoletzki, N., & Eyre-Walker, A. (2011). Estimation of the Neutrality Index. Molecular Biology and Evolution, 28(1), 63-70. https://doi.org/10.1093/molbev/msq249

## License

MIT License
