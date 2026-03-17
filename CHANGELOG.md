# Changelog

## [0.4.0] - 2026-03-17

### Added
- `mkado vcf` command: run MK tests directly from VCF + reference FASTA + GFF3 annotation
- Chunked parallel VCF processing for improved performance
- GFF3 parser with GTF format detection, malformed line handling, and hierarchy validation warnings
- `--verbose` flag for VCF mode to surface htslib, GFF parsing, and volcano plot diagnostics
- Example VCF dataset in `examples/example_vcf/`
- VCF input documentation (`docs/vcf-input.rst`)

### Fixed
- `cds_length` property called as method in GFF parser logging

## [0.3.0] - 2026-03-10

### Added
- Alternate genetic code support (`--code-table`) for 24 NCBI genetic code tables
- Name-based code selection (e.g. `--code-table vertebrate-mito`)
- `mkado codes` command to list available genetic code tables

## [0.2.0] - 2026-03-01

### Added
- Imputed MK test (Murga-Moreno et al. 2022) for correcting slightly deleterious mutations by imputation rather than discarding low-frequency polymorphisms
- `--imputed` flag for `test` and `batch` CLI commands
- Aggregated and per-gene batch modes for imputed MK test
- DFE decomposition (d, b, f fractions) when site counts are provided
- Documentation for the imputed MK test

## [0.1.2] - 2025-05-14

### Fixed
- Use delta method for `min_frequency` to correctly exclude singletons
- Remove title from volcano plots

### Added
- PyPI publish workflow (triggered by version tags)

## [0.1.1] - 2025-05-13

### Added
- Direction of Selection (DoS) statistic
- Tarone-Greenland alpha (α_TG) weighted multi-gene estimator
- `--no-singletons` option for automatic singleton exclusion
- Polarized MK test polymorphism polarization
- Volcano plot visualization for batch results
- Asymptotic alpha(x) plot visualization
- Benjamini-Hochberg adjusted p-values in batch output
- Sphinx documentation

## [0.1.0] - 2025-01-22

- Initial release

[0.4.0]: https://github.com/kr-colab/mkado/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/kr-colab/mkado/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/kr-colab/mkado/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/kr-colab/mkado/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/kr-colab/mkado/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/kr-colab/mkado/releases/tag/v0.1.0
