# Changelog

## [Unreleased]

### Fixed
- A codon pair where every mutational ordering passed through a stop codon
  was dropped from Dn/Ds/Pn/Ps with no signal, indistinguishable from an
  identical-codon pair. Under the vertebrate mitochondrial code this
  silently dropped four ordinary Lys/Trp pairs. The pair still drops from
  the counts, but each drop is now counted and surfaces as a warning, the
  same way skipped indels and multi-allelic sites already do (closes #90).
- Asymptotic MK frequency bin edges came from `np.linspace`, which does not
  always return the double nearest `k / num_bins` -- with the default 20
  bins, edge 3 came back as `0.15000000000000002` instead of `0.15`,
  misclassifying a polymorphism at exactly that frequency. Bin edges are
  now computed exactly (closes #84).
- Pooled-polymorphism frequency filtering was unreliable: `polarized_mk_test`'s
  pooled branch applied no filter at all, and `mk_test`'s pooled branch
  filtered on ingroup-only derived counts, missing outgroup-private
  singletons. Pooled mode now filters on the minor-allele frequency across
  both groups (closes #83).
- `alpha_tg_from_gene_data` reported a zero-width confidence interval for a
  single gene, or when every bootstrap replicate failed, indistinguishable
  from a genuinely tight CI. Both cases now report `ci_low`/`ci_high` as
  undefined (`NA` in pretty/TSV, `null` in JSON) instead, with `batch` and
  `vcf` emitting a warning (closes #79).
- `--imputed` combined with a polarize flag silently ran an unpolarized
  imputed test instead of erroring, unlike `--asymptotic` which already
  rejected the combination (with wording that had itself drifted
  between `test` and `batch`). `--plot-asymptotic`'s ignored-flag
  warning was missing its "(ignored)" suffix in `test`. `vcf`'s
  volcano-plot code caught a broader exception than `batch`'s identical
  call, swallowing real bugs. `--alpha-tg` combined with `--per-gene`
  silently ignored `--per-gene` with no warning. The `imputed_cutoff`
  and `is_aggregate` calculations were duplicated verbatim across
  `test`/`batch`/`vcf`. All six now go through shared helpers
  (closes #103).
- `test`, `batch`, and `vcf` each hand-wrote the same option-conflict
  checks and the copies had drifted: `vcf` dropped the explanatory
  clause from its messages, and checked `--alpha-tg`/`--imputed`
  against `--asymptotic` in a different order than `batch`, so the same
  flag combination could report a different conflict depending on the
  command. All three now share one validator (closes #74).
- Batch TSV chose its column layout from the first result's type, then
  silently dropped every row of a different type. A mixed-type batch
  now raises a clear error instead (closes #69).
- `--output -` (a bare dash for stdout) was rejected by the output-path
  check (closes #46).
- A failed VCF open leaked one file descriptor; the htslib capture pipe
  is now closed in every case (closes #45).
- `--no-singletons` did not filter singletons on the VCF path; the
  frequency check let the singleton itself through (closes #40).
- VCF divergence counting compared the outgroup codon to the reference
  instead of to the ingroup, miscounting sites where the ingroup was
  fixed for ALT (closes #44).
- `--freq-cutoffs` was not forwarded to per-gene asymptotic workers in
  `vcf` mode (closes #42).
- Batch TSV had no layout for imputed results and silently fell back to
  JSON (closes #43).
- Batch JSON silently dropped a gene when two results shared a name;
  duplicate names are now rejected up front (closes #70).
- CLI error messages had drifted: some exited without the "Error: "
  prefix. All CLI errors now route through one helper (closes #72).
- `batch` and `vcf` used different exit codes for an empty result set.
  Both now exit 0 when a run simply produces no results, and 1 only
  when genes errored (closes #75).
- `mkado vcf --gene X` with `-a`, `--imputed`, or `--alpha-tg` silently
  ran the standard MK test instead of the requested mode (closes #41).
- `--no-singletons` did not filter singletons on the FASTA
  `--polarize-match` or `--alpha-tg` paths, or wherever sample size did
  not match the raw sequence count. All FASTA paths now filter by
  derived-allele count, matching the VCF path (closes #55).
- Pooled major-codon selection weighted each population's frequency
  spectrum as a proportion, so a small population could outweigh a
  large one. Sequences are now counted directly (closes #81).
- Documented how a multi-position codon difference is resolved; the
  rule was previously unstated (closes #88).
- A multi-allelic codon collapsed two different derived bases into one
  mutation. Each base is now counted separately; polymorphism counts
  rise on affected genes (closes #87).
- VCF divergence rules now match the FASTA path: a codon touched by an
  ingroup polymorphism is excluded from divergence, and a heterozygous
  outgroup call is treated as missing rather than resolved. Divergence
  counts drop where this applies (closes #58).
- `--freq-cutoffs` had no effect on the FASTA per-gene asymptotic
  fitter; it always used the full frequency range (closes #61).
- Imputed batch results showed a placeholder adjusted p-value of 1.0
  instead of a real Benjamini-Hochberg correction; asymptotic per-gene
  batches showed a fabricated adjusted-p-value line despite having no
  p-value (closes #65).
- A failed curve fit in the asymptotic fitters fell back to a value
  outside the requested `--freq-cutoffs` window instead of the last
  point inside it.

## [0.5.0] - 2026-05-01

Round-1 revision feature release for the MKado applications note
(G3, MS G3-2026-406681). Adds the asymptotic-MK options reviewers
asked for, plus omega-decomposition reporting, performance work, and
documentation expansions; one small breaking change is the removal
of the Typer shell-completion options on the umbrella command.

### Added
- `--sfs-mode {at, above}` flag on `test`, `batch`, and `vcf` commands
  (closes #8). The `at` mode is the original per-bin SFS construction
  of Messer & Petrov (2013), in which `Pn(x)` and `Ps(x)` count
  polymorphic sites with derived allele frequency *equal to* x. The
  new `above` mode is the cumulative inclusive variant of Uricchio
  et al. (2019), in which `Pn(x)` and `Ps(x)` count polymorphic sites
  with derived allele frequency *at or above* x; the cumulative form
  has the same asymptote at x = 1 but is more robust to bin-count
  sparsity in large samples. Default remains `at` for backward
  compatibility.
- `--output / -O PATH` flag on `test`, `batch`, and `vcf` commands.
  Writes the formatted output to a named file rather than to stdout;
  progress bar and feedback messages remain on stderr. The
  `--format` flag continues to control human-readable / TSV / JSON.
- Nei-Gojobori site totals on every result dataclass (`MKResult`,
  `AsymptoticMKResult`, `PolarizedMKResult`, `ImputedMKResult`,
  `AlphaTGResult`) as `ln` and `ls` fields. Aggregated and per-gene
  paths populate them automatically from analyzed codons.
- `omega` (dN/dS, with Nei-Gojobori site weighting) on every result
  dataclass: `omega = (Dn/Ds) * (Ls/Ln)`. Returns `None` when `Ds = 0`
  or site counts are unavailable.
- `omega_a` and `omega_na` (adaptive / non-adaptive substitution rates)
  on result types whose alpha estimator is appropriate for the
  decomposition: `AsymptoticMKResult` (per-gene and aggregated),
  `ImputedMKResult`, and `AlphaTGResult`. Following Gossmann, Keightley
  & Eyre-Walker (2012, *Genome Biol Evol* 4(5):658-667),
  `omega_a = alpha * omega` and `omega_na = (1 - alpha) * omega`.
  `MKResult` and `PolarizedMKResult` deliberately omit these — the
  per-gene Smith-Eyre-Walker / polarized alpha is too noisy for a
  meaningful rate decomposition. See `docs/omega.rst` for the
  rationale.
- 95% confidence intervals on the omega decomposition where a sampling
  procedure is already available. `AsymptoticMKResult` exposes
  `omega_a_ci_low`/`omega_a_ci_high` and `omega_na_ci_low`/`omega_na_ci_high`
  by scaling the existing alpha CI (Dn/Ds/Ln/Ls are constants under the
  asymptotic Monte Carlo, so omega itself has no sampling distribution).
  `AlphaTGResult` recomputes omega/omega_a/omega_na per gene-resampled
  bootstrap replicate and reports `omega_ci_low/high`,
  `omega_a_ci_low/high`, `omega_na_ci_low/high` from the percentile
  distribution. `ImputedMKResult` does not (yet) bootstrap and reports
  point estimates only.
- `omega_decomposition()` helper in `mkado.analysis.statistics`.
- `AlignedPair.count_total_sites()` aggregating Nei-Gojobori `(Ln, Ls)`
  totals over analyzed codons.
- `Ln`, `Ls`, `omega`, `omega_a`, and `omega_na` columns in pretty,
  TSV, and JSON output for all result types.
- `--ci-method {monte-carlo,bootstrap}` CLI flag on `test`, `batch`, and
  `vcf` commands. The new `bootstrap` mode uses case-resampling of the
  pooled polymorphism list with replacement, refitting the curve per
  replicate (closes #10). Default remains `monte-carlo` to preserve
  existing behavior. Each result type now records the CI method via a
  `ci_method` field on `AsymptoticMKResult` (default `"monte-carlo"`),
  `AlphaTGResult` (`"bootstrap"`), and `ImputedMKResult`
  (`"bootstrap"` when set, `None` when no CI was computed).
- Bootstrap CI on imputed-MK alpha. `imputed_mk_test` and
  `imputed_mk_test_multi` now accept `n_bootstrap: int = 0` (preserves
  legacy behavior when 0) and `seed: int = 42`. CLI: enabled
  automatically when `--bootstrap > 0` (the existing flag drives the
  replicate count). New fields on `ImputedMKResult`:
  `alpha_ci_low/high`, `omega_a_ci_low/high`, `omega_na_ci_low/high`,
  `ci_method`. Following the asymptotic CI pattern, omega-decomposition
  CIs are derived analytically from the alpha CI scaled by omega.
- `_compute_ci_bootstrap()` helper in `mkado.analysis.asymptotic`
  alongside the existing `_compute_ci_monte_carlo()`.
- `_bootstrap_imputed_alpha()` helper in `mkado.analysis.imputed`.
- `examples/benchmarks/aggregated_asymptotic_runtime.py`
  (closes #11): reproducible runtime benchmark for
  `asymptotic_mk_test_aggregated` in the regime where the test is
  statistically reliable (~1,000 random genes per replicate, 100
  replicates per condition, full 2x2 grid of `sfs_mode` x
  `ci_method`). Polymorphism extraction is amortized via on-disk
  cache; figure caption and CSV header record the one-time
  extraction wall time so the cost composition is unambiguous.
  `--scaling` flag adds an end-to-end worker-count sweep
  (8/16/32/64). Documented at `examples/benchmarks/README.md`.
- `--workers N` exposed end-to-end on `test`, `batch`, and `vcf`
  (single CLI knob threads through to extraction and to the
  aggregated bootstrap-CI inner loop alike).
- "Working with Ortholog Alignments" tutorial section in
  `docs/tutorial.rst` documenting the OrthoMaM-style polymorphism-
  projection procedure used to assemble the human benchmark dataset
  (alignment matching, VCF extraction, haplotype assembly, alignment
  re-gapping, multi-transcript reconciliation).
- "Per-gene vs Aggregated Modes" section in `docs/batch-workflow.rst`
  explicitly contrasting `--aggregate` and `--per-gene`, naming the
  default for each test, and explaining why aggregation is the
  default for asymptotic MK.
- "Frequency-Threshold Correction (FWW)" subsection in
  `docs/alpha-tg.rst` showing how to combine `--alpha-tg` with
  `--min-freq X` to apply the Fay/Wyckoff/Wu correction within the
  Tarone-Greenland weighted estimate.
- "SFS Construction: At-x vs Above-x" section in `docs/asymptotic.rst`
  with the math and Uricchio 2019 citation.
- bioRxiv preprint badge in `README.md`.

### Performance
- `_bootstrap_imputed_alpha()` is now vectorized: pre-extracts numpy
  arrays of frequencies and N/S type once, then per replicate uses
  `rng.integers` indexing and `np.count_nonzero` instead of calling
  `_compute_imputed` (which does four Python-level `sum()` passes per
  call). ~63× speedup on the imputed-bootstrap path measured on a
  ~1000-polymorphism gene with 1000 replicates (1.4 s → 22 ms). CIs are
  byte-identical to the prior implementation under the same seed
  (closes #16).
- Per-gene `asymptotic_mk_test` bootstrap loop replaced by a call to
  the shared `_compute_ci_bootstrap()` helper (closes #17). Removes
  ~60 lines of duplicate code, vectorizes the per-replicate binning,
  and aligns failure-handling between per-gene and aggregated paths.
- Bootstrap CI inner loop parallelized via `ProcessPoolExecutor`
  with batched dispatch (closes #25). `asymptotic_mk_test_aggregated`
  with `--ci-method bootstrap` now scales across worker cores;
  threads were tried first but `scipy.optimize.curve_fit` runs the
  trust-region iteration in pure Python and holds the GIL.
- Analytic Jacobian for the asymptotic exponential and linear models
  (`_exponential_model_jac`, `_linear_model_jac`) replaces SciPy's
  finite-difference `approx_derivative` in `curve_fit`. Combined
  with vectorization of `aggregate_polymorphism_data` (using
  `np.searchsorted` + `np.bincount` over preallocated arrays in
  place of nested Python loops), this delivers ~28% end-to-end
  speedup on the aggregated-asymptotic-bootstrap workload.
  Numerically verified against the finite-difference path
  (rtol 1e-6 to 1e-9) in `tests/test_asymptotic.py::TestAnalyticJacobian`.

### Fixed
- `GeneticCode.translate` no longer uses `@lru_cache` on a bound method
  (closes #14). The decorator's cache key included `self`, pinning every
  `GeneticCode` instance for the lifetime of the class object — a real
  leak in tests, notebooks, and any process that constructs many
  ephemeral instances. The translation table is now eagerly precomputed
  in `__init__` as a per-instance dict (~640 bytes), mirroring the
  proven `_compute_codon_paths` pattern. Same fix applied to
  `count_synonymous_sites`, retiring its lazy-population dict for an
  eager precompute over the bounded 64-codon keyspace. Behavior is
  byte-identical for callers; a new regression test in
  `tests/test_codons.py::TestGeneticCodeMemory` uses `weakref` to verify
  instances are garbage-collected after `del` + `gc.collect()`.
- Standardized missing-value rendering on the literal `"NA"` (no slash)
  across pretty and TSV output (closes #13). The dataclass `__str__`
  methods previously rendered `"N/A"` while the TSV formatters used
  `"NA"` — a user piping pretty output through R/awk that parses `"NA"`
  would have silently missed `"N/A"` rows. JSON renders `null`,
  unchanged. New regression suite in `tests/test_output_format.py`
  guards both the convention and the absence of `"N/A"` in any output
  format.
- "Third outgroup" wording corrected to "second outgroup" in
  `docs/index.rst` and `README.md` (the polarized MK test uses two
  outgroups in total: the divergence outgroup and a polarization
  outgroup).
- `__version__` is now sourced from package metadata via
  `importlib.metadata.version("mkado")` instead of being hardcoded;
  `pyproject.toml` is the single source of truth, so the runtime
  version cannot drift from the released version. Module docstring
  corrected from the original "mikado" name.
- Per-gene `Pn`/`Ps` definitions tightened in result-dataclass
  docstrings to match the manuscript Methods (default zero
  threshold, `--min-freq`, `--no-singletons`, and the
  `--pool-polymorphisms` vs DnaSP convention spelled out).

### Changed
- Per-gene `asymptotic_mk_test` CI failure handling: replicates whose
  curve fit fails are now **dropped** rather than imputed with the last
  bin's alpha (the prior behavior biased the CI when fits failed
  systematically). When fewer than half of replicates succeed, the CI
  degenerates to the point estimate. Matches the aggregated path's
  behavior added in #10.

### Removed
- Typer's `--install-completion` and `--show-completion` umbrella
  options, by passing `add_completion=False` to the `Typer` app
  constructor. Reviewer 2 reported that completion was painfully
  slow on their system, which is consistent with `mkado`'s package
  import cost (scipy.stats + scipy.optimize together dominate
  start-up); every Tab press in completion mode pays this cost.
  Removing the feature is the cleaner answer than documenting it
  with a "may be slow" note. Users who genuinely want shell
  completion can still hand-roll it with their shell of choice.

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
