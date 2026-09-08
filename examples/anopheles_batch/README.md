# Anopheles Batch Example

This example demonstrates batch processing of MK tests using aligned coding sequences from two Anopheles mosquito species.

## Data

Four hundred randomly selected genes with aligned codon sequences (220 sequences each: 100 afun, 120 gamb). All genes have at least one observation in each cell of the MK contingency table (Dn, Ds, Pn, Ps >= 1).

Species:
- **afun**: Anopheles funestus (ingroup)
- **gamb**: Anopheles gambiae (outgroup)

## Example Commands

### Standard MK Test

```bash
# Batch MK test: afun (ingroup) vs gamb (outgroup)
mkado batch examples/anopheles_batch/ -i afun -o gamb

# TSV output for downstream analysis
mkado batch examples/anopheles_batch/ -i afun -o gamb -f tsv

# JSON output
mkado batch examples/anopheles_batch/ -i afun -o gamb -f json
```

### Asymptotic MK Test

```bash
# Aggregate asymptotic MK test (pools polymorphism across genes)
mkado batch examples/anopheles_batch/ -i afun -o gamb -a

# Per-gene asymptotic analysis
mkado batch examples/anopheles_batch/ -i afun -o gamb -a --per-gene

# Custom frequency bins
mkado batch examples/anopheles_batch/ -i afun -o gamb -a -b 5
```

### Single Gene Analysis

```bash
# Standard MK test on a single gene
mkado test examples/anopheles_batch/AGAP000078.fa -i afun -o gamb

# Asymptotic test on a single gene
mkado test examples/anopheles_batch/AGAP000078.fa -i afun -o gamb -a
```

## Expected Output

Running `mkado batch examples/anopheles_batch/ -i afun -o gamb` produces TSV output:

```
gene        Dn   Ds    Pn   Ps   p_value    p_value_adjusted  NI        alpha
AGAP000150  12   28    6    17   0.781067   0.949626          0.823529  0.176471
AGAP001182  31   44    19   32   0.712357   0.930941          0.842742  0.157258
...
```

Running the aggregate asymptotic MK test with 20 frequency bins fits an exponential model:

```
Dn     Ds     Pn    Ps     alpha_asymptotic  CI_low    CI_high   model        num_genes
18828  49857  4673  24340  0.690074          0.583694  0.718920  exponential  400
```

### Visualization Options

```bash
# Volcano plot (p-value vs alpha)
mkado batch examples/anopheles_batch/ -i afun -o gamb --volcano volcano.png

# Asymptotic alpha plot (aggregate)
mkado batch examples/anopheles_batch/ -i afun -o gamb -a --plot-asymptotic asymptotic.png
```

See the [batch workflow documentation](../../docs/batch-workflow.rst) for example plots.

## Notes

- Sequences are filtered by name pattern: `-i afun` matches all sequences containing "afun"
- The alignments are codon-aligned (in-frame); use `-r 1` (default) for reading frame
