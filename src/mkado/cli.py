"""Command-line interface for mkado."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, Optional

import typer
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.style import Style

from mkado import __version__
from mkado.analysis.asymptotic import (
    asymptotic_mk_test,
    asymptotic_mk_test_aggregated,
)
from mkado.analysis.mk_test import mk_test
from mkado.analysis.polarized import polarized_mk_test
from mkado.batch_workers import BatchTask, WorkerResult, _stop_blocked_warning, process_gene
from mkado.io.output import OutputFormat, format_batch_results, format_result
from scipy.stats import false_discovery_control

from mkado.analysis.mk_test import MKResult
from mkado.analysis.polarized import PolarizedMKResult

if TYPE_CHECKING:
    from mkado.analysis.alpha_tg import AlphaTGResult
    from mkado.io.output import BatchResult

# Console that writes to stderr (so progress doesn't mix with data output)
stderr_console = Console(stderr=True)


STDOUT_PATH = Path("-")


def cli_error(message: str) -> typer.Exit:
    """Report a user-facing error on stderr and return the exception to raise.

    Returning rather than raising keeps the ``raise`` at the call site, which is
    what typer.BadParameter already does in this module.
    """
    typer.echo(f"Error: {message}", err=True)
    return typer.Exit(1)


def resolve_code_table_or_exit(code_table: str) -> int:
    """Resolve a genetic code name or id, reporting an unknown one as a CLI error."""
    from mkado.data.genetic_codes import resolve_code_table

    try:
        return resolve_code_table(code_table)
    except ValueError as e:
        raise cli_error(str(e)) from e


def validate_path_not_flag(value: Path | None) -> Path | None:
    """Validate that a Path argument doesn't look like a flag.

    This catches common mistakes like: --option -a (where -a gets consumed as the path).
    A bare dash is a value, not an option, by POSIX and click convention, so it passes;
    the caller decides what it means.
    """
    if value is not None and value != STDOUT_PATH and str(value).startswith("-"):
        raise typer.BadParameter(
            f"'{value}' looks like a flag, not a file path. Check the order of your arguments."
        )
    return value


OutputOption = Annotated[
    Optional[Path],
    typer.Option(
        "--output",
        "-O",
        help="Write formatted results to this file (default: stdout). Use '-' for stdout.",
        callback=validate_path_not_flag,
    ),
]

CIMethodOption = Annotated[
    str,
    typer.Option(
        "--ci-method",
        help="CI method: 'monte-carlo' (default, parametric MVN sampling) or "
        "'bootstrap' (case-resampling). Meaningful for --asymptotic --aggregate; "
        "ignored elsewhere. Imputed CI is enabled automatically when --bootstrap > 0.",
    ),
]

SfsModeOption = Annotated[
    str,
    typer.Option(
        "--sfs-mode",
        help="Asymptotic-MK SFS construction: 'at' (default, Messer & Petrov 2013, "
        "per-bin counts) or 'above' (Uricchio et al. 2019, inclusive right-tail "
        "cumulative counts). Both modes share the asymptote at x=1 but 'above' "
        "is more sample-size-stable. Ignored when --asymptotic is not set.",
    ),
]


def resolve_output_format(output_format: str) -> OutputFormat:
    """Convert the --format string to its enum member, reporting an unknown one."""
    try:
        return OutputFormat(output_format)
    except ValueError as e:
        raise cli_error(f"Invalid format '{output_format}'.") from e


def validate_ci_method(ci_method: str) -> None:
    """Reject any ``--ci-method`` value other than the two supported strings."""
    if ci_method not in ("monte-carlo", "bootstrap"):
        raise cli_error(f"Invalid --ci-method '{ci_method}'. Use 'monte-carlo' or 'bootstrap'.")


def validate_sfs_mode(sfs_mode: str) -> None:
    """Reject any ``--sfs-mode`` value other than 'at' or 'above'."""
    if sfs_mode not in ("at", "above"):
        raise cli_error(f"Invalid --sfs-mode '{sfs_mode}'. Use 'at' or 'above'.")


def validate_option_compatibility(
    *,
    use_asymptotic: bool,
    use_imputed: bool,
    no_singletons: bool,
    min_freq: float,
    alpha_tg: bool = False,
) -> None:
    """Reject combinations of frequency-filter and alpha-method flags that conflict.

    Shared by ``test``, ``batch``, and ``vcf`` so all three raise the same message in
    the same precedence order. ``alpha_tg`` defaults to False for callers without an
    ``--alpha-tg`` flag; ``test`` has none and passes False explicitly.
    """
    freq_cutoffs_note = "The asymptotic test uses --freq-cutoffs for frequency filtering."
    if use_asymptotic and min_freq > 0.0:
        raise cli_error(f"--min-freq cannot be used with --asymptotic. {freq_cutoffs_note}")
    if alpha_tg and use_asymptotic:
        raise cli_error(
            "--alpha-tg and --asymptotic are mutually exclusive. "
            "Choose one method for estimating alpha."
        )
    if use_asymptotic and no_singletons:
        raise cli_error(f"--no-singletons cannot be used with --asymptotic. {freq_cutoffs_note}")
    if use_imputed and use_asymptotic:
        raise cli_error("--imputed and --asymptotic are mutually exclusive.")
    if use_imputed and alpha_tg:
        raise cli_error("--imputed and --alpha-tg are mutually exclusive.")
    if use_imputed and no_singletons:
        raise cli_error(
            "--no-singletons cannot be used with --imputed. "
            "The imputed test needs low-frequency variants."
        )
    if no_singletons and min_freq > 0.0:
        raise cli_error(
            "--no-singletons and --min-freq cannot be used together. "
            "--no-singletons filters by allele count instead."
        )


def validate_polarize_compatibility(
    *, use_asymptotic: bool, use_imputed: bool, polarize: bool, polarize_flag: str
) -> None:
    """Reject a second-outgroup flag combined with a test mode that can't use it.

    ``polarize_flag`` names the flag actually in play, since it differs by command
    and mode: ``--polarize-match``, ``-p/--polarize``, or ``--polarize-pattern``.
    """
    if use_asymptotic and polarize:
        raise cli_error(f"--asymptotic and {polarize_flag} are mutually exclusive.")
    if use_imputed and polarize:
        raise cli_error(f"--imputed and {polarize_flag} are mutually exclusive.")


def resolve_imputed_cutoff(*, use_imputed: bool, min_freq: float) -> float:
    """Resolve the imputed test's derived-allele-frequency cutoff from --min-freq."""
    return min_freq if (use_imputed and min_freq > 0.0) else 0.15


def resolve_is_aggregate(
    *, use_asymptotic: bool, use_imputed: bool, alpha_tg: bool, aggregate: bool
) -> bool:
    """alpha_TG is always aggregate; asymptotic/imputed follow --aggregate/--per-gene."""
    return (use_asymptotic and aggregate) or alpha_tg or (use_imputed and aggregate)


def warn_if_plot_asymptotic_ignored(*, plot_asymptotic: Path | None, use_asymptotic: bool) -> None:
    """Warn that --plot-asymptotic has no effect without --asymptotic."""
    if plot_asymptotic and not use_asymptotic:
        typer.echo("Warning: --plot-asymptotic requires --asymptotic/-a flag (ignored)", err=True)


def warn_if_alpha_tg_ignores_per_gene(*, alpha_tg: bool, aggregate: bool) -> None:
    """Warn that --per-gene has no effect on --alpha-tg, which is always aggregate."""
    if alpha_tg and not aggregate:
        typer.echo(
            "Warning: --alpha-tg always aggregates across genes; --per-gene is ignored",
            err=True,
        )


def warn_if_alpha_tg_ci_undefined(result: AlphaTGResult) -> None:
    """Warn that the alpha_TG bootstrap CI is undefined with fewer than two genes."""
    if result.num_genes < 2:
        typer.echo(
            "Warning: alpha_TG confidence interval requires at least 2 genes "
            f"(got {result.num_genes}); CI fields reported as NA",
            err=True,
        )


def parse_frequency_cutoffs(freq_cutoffs: str) -> tuple[float, float]:
    """Parse a ``--freq-cutoffs 'low,high'`` string into a (low, high) tuple."""
    try:
        low, high = freq_cutoffs.split(",")
        return (float(low), float(high))
    except ValueError:
        raise cli_error(f"Invalid frequency cutoffs '{freq_cutoffs}'. Use 'low,high'")


def resolve_ci_replicates(bootstrap: int, ci_method: str) -> int:
    """Map ``--bootstrap`` to ``ci_replicates`` for the chosen CI method.

    MC samples params from a covariance matrix (cheap, ~10000×); bootstrap
    refits the curve per replicate (expensive, ~100×). The 100× factor matches
    the cost ratio between the two paths.
    """
    return bootstrap * 100 if ci_method == "monte-carlo" else bootstrap


def write_output(content: str, output_path: Path | None) -> None:
    """Write formatted result content to a file, or to stdout if no path given.

    A path of ``None`` or ``Path("-")`` means stdout (POSIX convention).
    """
    if output_path is None or output_path == STDOUT_PATH:
        typer.echo(content)
        return
    text = content if content.endswith("\n") else content + "\n"
    output_path.write_text(text)
    typer.echo(f"Results saved to {output_path}", err=True)


def _collect_gene_data(
    worker_results: list[WorkerResult],
    mode_label: str,
) -> list:
    """Unpack the per-gene results and report how many genes contributed.

    ``run_parallel_batch`` keeps only results that carry data, so every entry
    here counts.

    Args:
        worker_results: List of WorkerResult from parallel processing.
        mode_label: Label for the stderr message (e.g., "aggregated asymptotic").

    Returns:
        List of result objects.
    """
    gene_data = [r.result for r in worker_results]
    typer.echo(f"Using {len(gene_data)} genes for {mode_label}", err=True)
    return gene_data


def warn_duplicate_stems(files: list[Path]) -> None:
    """Warn when input files share a stem, because the stem becomes the gene name."""
    by_stem: dict[str, list[Path]] = {}
    for path in files:
        by_stem.setdefault(path.stem, []).append(path)
    for stem, paths in sorted(by_stem.items()):
        if len(paths) > 1:
            names = ", ".join(str(path) for path in paths)
            typer.echo(
                f"Warning: {len(paths)} files share the gene name '{stem}': {names}", err=True
            )


def write_batch_output(
    results: list[tuple[str, BatchResult]],
    fmt: OutputFormat,
    adjusted_pvalues: list[float] | None,
    output: Path | None,
) -> None:
    """Format batch results and write them, reporting a rejected result set as an error."""
    try:
        content = format_batch_results(results, fmt, adjusted_pvalues)
    except ValueError as e:
        raise cli_error(str(e)) from e
    write_output(content, output)


def save_volcano_plot(results: list[tuple[str, BatchResult]], path: Path) -> None:
    """Write a volcano plot for batch/vcf results, or report why it couldn't be made."""
    from mkado.io.plotting import create_volcano_plot

    try:
        create_volcano_plot(results, path)
        typer.echo(f"Volcano plot saved to {path}", err=True)
    except ValueError as e:
        typer.echo(f"Could not generate volcano plot: {e}", err=True)


def compute_adjusted_pvalues(
    results: list[tuple[str, BatchResult]],
) -> list[float] | None:
    """Compute Benjamini-Hochberg adjusted p-values for batch results.

    Args:
        results: List of (name, result) tuples from MK tests

    Returns:
        Adjusted p-values in the same order as ``results``, or ``None`` when
        the result type carries no p-value (currently only
        ``AsymptoticMKResult``) or ``results`` is empty.
    """
    from mkado.analysis.asymptotic import AsymptoticMKResult
    from mkado.analysis.imputed import ImputedMKResult

    p_values = []
    for _, result in results:
        if isinstance(result, MKResult):
            p_values.append(result.p_value)
        elif isinstance(result, PolarizedMKResult):
            p_values.append(result.p_value_ingroup)
        elif isinstance(result, ImputedMKResult):
            p_values.append(result.p_value)
        elif isinstance(result, AsymptoticMKResult):
            # Every call site batches one result type per invocation, so
            # dropping these entries still leaves p_values index-aligned
            # with `results` (all-or-nothing: either every entry is skipped
            # here and the empty-p_values check below returns None, or none
            # is).
            continue
        else:
            raise TypeError(f"Unknown result type: {type(result)}")

    if not p_values:
        return None

    return list(false_discovery_control(p_values, method="bh"))


class RainbowBarColumn(BarColumn):
    """A progress bar that cycles through rainbow colors."""

    RAINBOW_COLORS = [
        "#FF0000",  # Red
        "#FF7F00",  # Orange
        "#FFFF00",  # Yellow
        "#00FF00",  # Green
        "#0000FF",  # Blue
        "#4B0082",  # Indigo
        "#9400D3",  # Violet
    ]

    def __init__(self) -> None:
        super().__init__(bar_width=40)
        self._color_index = 0

    def render(self, task):  # type: ignore[no-untyped-def]
        """Render the bar with rainbow colors."""
        if task.total:
            progress = task.completed / task.total
            color_idx = int(progress * len(self.RAINBOW_COLORS) * 3) % len(self.RAINBOW_COLORS)
        else:
            color_idx = self._color_index
            self._color_index = (self._color_index + 1) % len(self.RAINBOW_COLORS)

        self.complete_style = Style(color=self.RAINBOW_COLORS[color_idx])
        self.finished_style = Style(color="#9400D3")
        return super().render(task)


def create_rainbow_progress() -> Progress:
    """Create a rainbow-colored progress bar that writes to stderr."""
    return Progress(
        SpinnerColumn(style="bold magenta"),
        TextColumn("[bold blue]{task.description}"),
        RainbowBarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=stderr_console,
    )


def get_worker_count(requested: int, num_tasks: int) -> int:
    """Determine the optimal number of workers."""
    if requested == 1 or num_tasks < 10:
        return 1
    cpu_count = os.cpu_count() or 4
    if requested > 0:
        return min(requested, cpu_count)
    return max(1, min(cpu_count - 1, num_tasks))


def report_no_results(had_error: bool) -> None:
    """Report an empty result set.

    Producing no results is a normal outcome, so it succeeds. It is a failure only
    when the set is empty because genes errored, which the per-gene messages detail.
    """
    if had_error:
        raise cli_error("No results to display; every gene failed")
    typer.echo("No results to display", err=True)


def run_parallel_batch(
    tasks: list[BatchTask],
    num_workers: int,
    description: str,
) -> tuple[list[WorkerResult], list[str], bool]:
    """Run batch processing with ProcessPoolExecutor.

    Returns the successful results, the messages to print, and whether any task errored.
    """
    results: list[WorkerResult] = []
    warnings: list[str] = []
    had_error = False

    if num_workers == 1:
        with create_rainbow_progress() as progress:
            task_id = progress.add_task(description, total=len(tasks))
            for task in tasks:
                worker_result = process_gene(task)
                if worker_result.error:
                    warnings.append(worker_result.error)
                    had_error = True
                elif worker_result.warning:
                    warnings.append(f"Warning: {worker_result.warning}")
                elif worker_result.result is not None:
                    results.append(worker_result)
                progress.advance(task_id)
    else:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_gene, task): task for task in tasks}
            with create_rainbow_progress() as progress:
                task_id = progress.add_task(description, total=len(tasks))
                for future in as_completed(futures):
                    try:
                        worker_result = future.result()
                        if worker_result.error:
                            warnings.append(worker_result.error)
                            had_error = True
                        elif worker_result.warning:
                            warnings.append(f"Warning: {worker_result.warning}")
                        elif worker_result.result is not None:
                            results.append(worker_result)
                    except Exception as e:
                        task = futures[future]
                        warnings.append(f"Error processing {task.file_path.name}: {e}")
                        had_error = True
                    progress.advance(task_id)

    return results, warnings, had_error


def find_partner_file(ingroup_file: Path, candidates: list[Path]) -> Path | None:
    """Pick the candidate that belongs to the same gene as an ingroup file.

    The gene name is the ingroup stem without its ``_ingroup`` or ``_in``
    marker. A candidate belongs to the gene when its stem is that name
    followed by ``_``: a bare prefix would let ``g1`` claim ``g10``'s file.
    The first qualifying candidate wins, so the caller passes them sorted and
    the filesystem's listing order does not decide between two that qualify.
    """
    stem = ingroup_file.stem
    gene = stem.removesuffix("_ingroup") if stem.endswith("_ingroup") else stem.removesuffix("_in")
    for candidate in candidates:
        if candidate.stem.startswith(gene + "_"):
            return candidate
    return None


def find_alignment_files(input_dir: Path) -> list[Path]:
    """Auto-detect alignment files in a directory.

    Tries common FASTA extensions in order: *.fa, *.fasta, *.fna
    Returns files matching the first pattern that finds results.
    """
    for pattern in ["*.fa", "*.fasta", "*.fna"]:
        files = sorted(input_dir.glob(pattern))
        if files:
            return files
    return []


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        typer.echo(f"mkado {__version__}")
        raise typer.Exit()


app = typer.Typer(
    name="mkado",
    help="MKado 御門: McDonald-Kreitman test toolkit.\n\n"
    "A modern Python implementation for detecting selection using "
    "the McDonald-Kreitman test and related methods.",
    no_args_is_help=True,
    add_completion=False,
)


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option("--version", "-v", callback=version_callback, is_eager=True),
    ] = False,
) -> None:
    """mkado: McDonald-Kreitman test toolkit."""
    pass


@app.command()
def test(
    fasta: Annotated[
        Path,
        typer.Argument(help="FASTA file (combined alignment, or ingroup sequences)"),
    ],
    outgroup_file: Annotated[
        Optional[Path],
        typer.Argument(help="Outgroup FASTA file (only for separate files mode)"),
    ] = None,
    # === Sequence filtering (combined file mode) ===
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--ingroup-match",
            "-i",
            help="Ingroup sequence name pattern (enables combined file mode)",
        ),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--outgroup-match",
            "-o",
            help="Outgroup sequence name pattern (required with -i)",
        ),
    ] = None,
    polarize_match: Annotated[
        Optional[str],
        typer.Option(
            "--polarize-match",
            help="Second outgroup pattern for polarized test (combined mode)",
        ),
    ] = None,
    # === Polarization (separate files mode) ===
    polarize_file: Annotated[
        Optional[Path],
        typer.Option(
            "--polarize",
            "-p",
            help="Second outgroup file for polarized test (separate files mode)",
        ),
    ] = None,
    # === Analysis type ===
    use_asymptotic: Annotated[
        bool,
        typer.Option(
            "--asymptotic",
            "-a",
            help="Use asymptotic MK test (accounts for slightly deleterious mutations)",
        ),
    ] = False,
    use_imputed: Annotated[
        bool,
        typer.Option(
            "--imputed",
            help="Use imputed MK test (Murga-Moreno et al. 2022). "
            "Uses --min-freq as DAF cutoff (default 0.15 if not set).",
        ),
    ] = False,
    # === Asymptotic options ===
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Frequency bins for asymptotic test"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option("--bootstrap", help="Bootstrap replicates for CI (asymptotic)"),
    ] = 100,
    ci_method: CIMethodOption = "monte-carlo",
    sfs_mode: SfsModeOption = "at",
    freq_cutoffs: Annotated[
        str,
        typer.Option("--freq-cutoffs", help="Frequency range 'low,high' (asymptotic)"),
    ] = "0.1,0.9",
    workers: Annotated[
        int,
        typer.Option(
            "--workers",
            "-w",
            min=0,
            help=(
                "Parallel workers for the bootstrap CI inner loop "
                "(0=auto, 1=sequential). Only meaningful for --asymptotic "
                "with --ci-method bootstrap; ignored otherwise."
            ),
        ),
    ] = 1,
    # === Common options ===
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format: pretty, tsv, json"),
    ] = "pretty",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1-3)"),
    ] = 1,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option("--pool-polymorphisms", help="Pool polymorphisms from both populations"),
    ] = False,
    min_freq: Annotated[
        float,
        typer.Option("--min-freq", help="Min derived allele frequency (0.0-1.0)"),
    ] = 0.0,
    no_singletons: Annotated[
        bool,
        typer.Option(
            "--no-singletons",
            help="Exclude sites where the derived allele appears once",
        ),
    ] = False,
    code_table: Annotated[
        str,
        typer.Option(
            "--code-table",
            help="Genetic code: name (e.g. vertebrate-mito) or NCBI table ID. "
            "Run 'mkado codes' to list options.",
        ),
    ] = "standard",
    plot_asymptotic: Annotated[
        Optional[Path],
        typer.Option(
            "--plot-asymptotic",
            help="Generate alpha(x) plot for asymptotic test (PNG, PDF, or SVG)",
            callback=validate_path_not_flag,
        ),
    ] = None,
    output: OutputOption = None,
) -> None:
    """Run McDonald-Kreitman test on a single alignment.

    TWO MODES OF OPERATION:

    1. COMBINED FILE MODE (recommended):
       Use -i and -o to filter sequences by name pattern from a single alignment.

       mkado test alignment.fa -i "speciesA" -o "speciesB"
       mkado test alignment.fa -i "gamb" -o "002019" --asymptotic

    2. SEPARATE FILES MODE:
       Provide two FASTA files (ingroup and outgroup).

       mkado test ingroup.fa outgroup.fa
       mkado test ingroup.fa outgroup.fa -p outgroup2.fa  # polarized

    ANALYSIS TYPES:

    - Standard MK test (default): Classic McDonald-Kreitman test
    - Asymptotic MK test (-a): Accounts for slightly deleterious mutations
      by examining alpha across the frequency spectrum (Messer & Petrov 2013)

    EXAMPLES:

        mkado test alignment.fa -i "dmel" -o "dsim"
        mkado test alignment.fa -i "dmel" -o "dsim" -a -b 20
        mkado test alignment.fa -i "dmel" -o "dsim" -a --freq-cutoffs 0.2,0.8
        mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"
        mkado test ingroup.fa outgroup.fa
        mkado test ingroup.fa outgroup.fa -a
    """
    from mkado.core.sequences import SequenceSet

    fmt = resolve_output_format(output_format)

    validate_ci_method(ci_method)
    validate_sfs_mode(sfs_mode)

    frequency_cutoffs = parse_frequency_cutoffs(freq_cutoffs)

    validate_option_compatibility(
        use_asymptotic=use_asymptotic,
        use_imputed=use_imputed,
        no_singletons=no_singletons,
        min_freq=min_freq,
        alpha_tg=False,
    )

    imputed_cutoff = resolve_imputed_cutoff(use_imputed=use_imputed, min_freq=min_freq)

    # Build genetic code
    from mkado.core.codons import GeneticCode

    code_table_id = resolve_code_table_or_exit(code_table)

    genetic_code = GeneticCode(table_id=code_table_id) if code_table_id != 1 else None

    # Determine mode
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            raise cli_error(
                "Combined mode requires both -i/--ingroup-match and -o/--outgroup-match"
            )

        if outgroup_file is not None:
            raise cli_error("Don't provide outgroup file when using -i/-o (combined mode)")

        if polarize_file is not None:
            raise cli_error("Use --polarize-match instead of -p in combined mode")

        validate_polarize_compatibility(
            use_asymptotic=use_asymptotic,
            use_imputed=use_imputed,
            polarize=bool(polarize_match),
            polarize_flag="--polarize-match",
        )

        # Load and filter sequences
        all_seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            raise cli_error(f"No sequences match ingroup pattern '{ingroup_match}'")
        if len(outgroup_seqs) == 0:
            raise cli_error(f"No sequences match outgroup pattern '{outgroup_match}'")

        typer.echo(
            f"Found {len(ingroup_seqs)} ingroup, {len(outgroup_seqs)} outgroup sequences",
            err=True,
        )

        if no_singletons:
            typer.echo("Excluding singletons", err=True)

        # Run appropriate test
        if use_asymptotic:
            result = asymptotic_mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                num_bins=bins,
                bootstrap_replicates=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
                genetic_code=genetic_code,
                sfs_mode=sfs_mode,
                frequency_cutoffs=frequency_cutoffs,
                workers=get_worker_count(workers, max(bootstrap, 1)),
            )
        elif use_imputed:
            from mkado.analysis.asymptotic import extract_polymorphism_data
            from mkado.analysis.imputed import imputed_mk_test

            poly_data = extract_polymorphism_data(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                genetic_code=genetic_code,
            )
            result = imputed_mk_test(poly_data, cutoff=imputed_cutoff, n_bootstrap=bootstrap)
        elif polarize_match:
            outgroup2_seqs = all_seqs.filter_by_name(polarize_match)
            if len(outgroup2_seqs) == 0:
                raise cli_error(f"No sequences match polarize pattern '{polarize_match}'")
            typer.echo(f"Polarizing with {len(outgroup2_seqs)} outgroup2 sequences", err=True)
            result = polarized_mk_test(
                ingroup=ingroup_seqs,
                outgroup1=outgroup_seqs,
                outgroup2=outgroup2_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                no_singletons=no_singletons,
                genetic_code=genetic_code,
            )
        else:
            result = mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                no_singletons=no_singletons,
                genetic_code=genetic_code,
            )

    else:
        # Separate files mode
        if outgroup_file is None:
            raise cli_error("Provide outgroup file, or use -i/-o for combined file mode")

        if polarize_match is not None:
            raise cli_error("Use -p/--polarize instead of --polarize-match in separate files mode")

        validate_polarize_compatibility(
            use_asymptotic=use_asymptotic,
            use_imputed=use_imputed,
            polarize=bool(polarize_file),
            polarize_flag="-p/--polarize",
        )

        if no_singletons:
            typer.echo("Excluding singletons", err=True)

        # Run appropriate test
        if use_asymptotic:
            result = asymptotic_mk_test(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                num_bins=bins,
                bootstrap_replicates=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
                genetic_code=genetic_code,
                sfs_mode=sfs_mode,
                frequency_cutoffs=frequency_cutoffs,
                workers=get_worker_count(workers, max(bootstrap, 1)),
            )
        elif use_imputed:
            from mkado.analysis.asymptotic import extract_polymorphism_data
            from mkado.analysis.imputed import imputed_mk_test

            poly_data = extract_polymorphism_data(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                genetic_code=genetic_code,
            )
            result = imputed_mk_test(poly_data, cutoff=imputed_cutoff, n_bootstrap=bootstrap)
        elif polarize_file:
            result = polarized_mk_test(
                ingroup=fasta,
                outgroup1=outgroup_file,
                outgroup2=polarize_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                no_singletons=no_singletons,
                genetic_code=genetic_code,
            )
        else:
            result = mk_test(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                no_singletons=no_singletons,
                genetic_code=genetic_code,
            )

    stop_blocked_warning = _stop_blocked_warning(fasta.stem, result)
    if stop_blocked_warning:
        typer.echo(f"Warning: {stop_blocked_warning}", err=True)

    write_output(format_result(result, fmt), output)

    # Generate asymptotic plot if requested
    if plot_asymptotic and use_asymptotic:
        from mkado.analysis.asymptotic import AsymptoticMKResult
        from mkado.io.plotting import create_asymptotic_plot

        if isinstance(result, AsymptoticMKResult):
            try:
                create_asymptotic_plot(result, plot_asymptotic)
                typer.echo(f"Asymptotic plot saved to {plot_asymptotic}", err=True)
            except ValueError as e:
                typer.echo(f"Could not generate plot: {e}", err=True)
    warn_if_plot_asymptotic_ignored(plot_asymptotic=plot_asymptotic, use_asymptotic=use_asymptotic)


@app.command()
def batch(
    input_dir: Annotated[
        Path,
        typer.Argument(help="Directory containing FASTA alignment files"),
    ],
    # === Sequence filtering (combined file mode) ===
    ingroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--ingroup-match",
            "-i",
            help="Ingroup sequence name pattern (enables combined file mode)",
        ),
    ] = None,
    outgroup_match: Annotated[
        Optional[str],
        typer.Option(
            "--outgroup-match",
            "-o",
            help="Outgroup sequence name pattern (required with -i)",
        ),
    ] = None,
    polarize_match: Annotated[
        Optional[str],
        typer.Option(
            "--polarize-match",
            help="Second outgroup pattern for polarized test",
        ),
    ] = None,
    # === File patterns ===
    file_pattern: Annotated[
        Optional[str],
        typer.Option(
            "--pattern",
            help="File glob pattern (default: auto-detect *.fa, *.fasta, *.fna)",
        ),
    ] = None,
    # === Separate files mode options ===
    ingroup_pattern: Annotated[
        str,
        typer.Option(help="Ingroup file pattern (separate files mode)"),
    ] = "*_ingroup.fa",
    outgroup_pattern: Annotated[
        str,
        typer.Option(help="Outgroup file pattern (separate files mode)"),
    ] = "*_outgroup.fa",
    polarize_pattern: Annotated[
        Optional[str],
        typer.Option(
            "--polarize-pattern",
            help="Second outgroup file pattern (separate files polarized mode)",
        ),
    ] = None,
    # === Analysis type ===
    use_asymptotic: Annotated[
        bool,
        typer.Option(
            "--asymptotic",
            "-a",
            help="Use asymptotic MK test",
        ),
    ] = False,
    aggregate: Annotated[
        bool,
        typer.Option(
            "--aggregate/--per-gene",
            help="Aggregate across genes (default) or per-gene results (asymptotic)",
        ),
    ] = True,
    alpha_tg: Annotated[
        bool,
        typer.Option(
            "--alpha-tg",
            help="Compute α_TG (Stoletzki & Eyre-Walker 2011 weighted estimator)",
        ),
    ] = False,
    use_imputed: Annotated[
        bool,
        typer.Option(
            "--imputed",
            help="Use imputed MK test (Murga-Moreno et al. 2022). "
            "Uses --min-freq as DAF cutoff (default 0.15 if not set).",
        ),
    ] = False,
    # === Asymptotic options ===
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Frequency bins (asymptotic)"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option("--bootstrap", help="Bootstrap replicates (asymptotic)"),
    ] = 100,
    ci_method: CIMethodOption = "monte-carlo",
    sfs_mode: SfsModeOption = "at",
    freq_cutoffs: Annotated[
        str,
        typer.Option("--freq-cutoffs", help="Frequency range 'low,high' (asymptotic)"),
    ] = "0.1,0.9",
    # === Common options ===
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format: pretty, tsv, json"),
    ] = "tsv",
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1-3)"),
    ] = 1,
    pool_polymorphisms: Annotated[
        bool,
        typer.Option("--pool-polymorphisms", help="Pool polymorphisms from both populations"),
    ] = False,
    min_freq: Annotated[
        float,
        typer.Option("--min-freq", help="Min derived allele frequency"),
    ] = 0.0,
    no_singletons: Annotated[
        bool,
        typer.Option(
            "--no-singletons",
            help="Exclude sites where the derived allele appears once",
        ),
    ] = False,
    code_table: Annotated[
        str,
        typer.Option(
            "--code-table",
            help="Genetic code: name (e.g. vertebrate-mito) or NCBI table ID. "
            "Run 'mkado codes' to list options.",
        ),
    ] = "standard",
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", min=0, help="Parallel workers (0=auto, 1=sequential)"),
    ] = 0,
    volcano: Annotated[
        Optional[Path],
        typer.Option(
            "--volcano",
            help="Generate volcano plot and save to specified path (PNG, PDF, or SVG)",
            callback=validate_path_not_flag,
        ),
    ] = None,
    plot_asymptotic: Annotated[
        Optional[Path],
        typer.Option(
            "--plot-asymptotic",
            help="Generate alpha(x) plot for aggregated asymptotic test (PNG, PDF, or SVG)",
            callback=validate_path_not_flag,
        ),
    ] = None,
    output: OutputOption = None,
) -> None:
    """Run MK test on multiple alignment files.

    TWO MODES OF OPERATION:

    1. COMBINED FILE MODE (recommended):
       Each file contains sequences from multiple species.
       Use -i and -o to filter by sequence name pattern.
       Mode is auto-detected when -i is provided.

       mkado batch alignments/ -i "speciesA" -o "speciesB"

    2. SEPARATE FILES MODE:
       Pairs of ingroup/outgroup files (e.g., gene1_ingroup.fa, gene1_outgroup.fa).
       Used when -i is NOT provided.

       mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"

    FILE DETECTION:
       Files are auto-detected by trying *.fa, *.fasta, *.fna in order.
       Override with --pattern "*.your_extension"

    ANALYSIS TYPES:

    - Standard MK test (default): Per-gene Dn, Ds, Pn, Ps
    - Asymptotic MK (-a): Genome-wide alpha with frequency correction
      - --aggregate (default): Pool data across genes, fit single curve
      - --per-gene: Fit curve for each gene separately

    EXAMPLES:

        mkado batch alignments/ -i "dmel" -o "dsim"
        mkado batch alignments/ -i "dmel" -o "dsim" -a
        mkado batch alignments/ -i "dmel" -o "dsim" -a --per-gene
        mkado batch alignments/ -i "dmel" -o "dsim" -w 8
        mkado batch genes/ --ingroup-pattern "*_in.fa" --outgroup-pattern "*_out.fa"
    """
    validate_ci_method(ci_method)
    validate_sfs_mode(sfs_mode)

    fmt = resolve_output_format(output_format)

    # Resolve genetic code table
    code_table_id = resolve_code_table_or_exit(code_table)

    validate_option_compatibility(
        use_asymptotic=use_asymptotic,
        use_imputed=use_imputed,
        no_singletons=no_singletons,
        min_freq=min_freq,
        alpha_tg=alpha_tg,
    )

    imputed_cutoff = resolve_imputed_cutoff(use_imputed=use_imputed, min_freq=min_freq)

    frequency_cutoffs = parse_frequency_cutoffs(freq_cutoffs)

    is_aggregate = resolve_is_aggregate(
        use_asymptotic=use_asymptotic,
        use_imputed=use_imputed,
        alpha_tg=alpha_tg,
        aggregate=aggregate,
    )
    warn_if_alpha_tg_ignores_per_gene(alpha_tg=alpha_tg, aggregate=aggregate)

    # Auto-detect mode based on -i flag
    combined_mode = ingroup_match is not None

    if combined_mode:
        # === COMBINED FILE MODE ===
        if not outgroup_match:
            raise cli_error("-o/--outgroup-match required with -i/--ingroup-match")

        validate_polarize_compatibility(
            use_asymptotic=use_asymptotic,
            use_imputed=use_imputed,
            polarize=bool(polarize_match),
            polarize_flag="--polarize-match",
        )

        # Find alignment files
        if file_pattern:
            alignment_files = sorted(input_dir.glob(file_pattern))
        else:
            alignment_files = find_alignment_files(input_dir)

        if not alignment_files:
            pattern_msg = file_pattern or "*.fa, *.fasta, *.fna"
            raise cli_error(f"No files found in {input_dir} (tried: {pattern_msg})")

        typer.echo(f"Found {len(alignment_files)} alignment files", err=True)
        warn_duplicate_stems(alignment_files)

        # Determine workers
        num_workers = get_worker_count(workers, len(alignment_files))

        # Build tasks
        tasks = [
            BatchTask(
                file_path=f,
                ingroup_match=ingroup_match,
                outgroup_match=outgroup_match,
                polarize_match=polarize_match,
                reading_frame=reading_frame,
                use_asymptotic=use_asymptotic,
                use_imputed=use_imputed and not aggregate,
                imputed_cutoff=imputed_cutoff,
                bins=bins,
                bootstrap=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
                min_freq=min_freq,
                no_singletons=no_singletons,
                extract_only=is_aggregate,
                code_table=code_table_id,
                ci_method=ci_method,
                sfs_mode=sfs_mode,
                frequency_cutoffs=frequency_cutoffs,
            )
            for f in alignment_files
        ]

        # Alpha TG mode
        if alpha_tg:
            from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings, _ = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = _collect_gene_data(worker_results, "alpha-TG")

            if gene_data_list:
                result = alpha_tg_from_gene_data(
                    gene_data=gene_data_list,
                    bootstrap_replicates=bootstrap,
                )
                warn_if_alpha_tg_ci_undefined(result)
                write_output(format_result(result, fmt), output)
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Aggregated asymptotic mode
        if use_asymptotic and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings, _ = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = _collect_gene_data(
                worker_results, "aggregated asymptotic"
            )

            if gene_data_list:
                ci_replicates = resolve_ci_replicates(bootstrap, ci_method)
                result = asymptotic_mk_test_aggregated(
                    gene_data=gene_data_list,
                    num_bins=bins,
                    ci_replicates=ci_replicates,
                    frequency_cutoffs=frequency_cutoffs,
                    ci_method=ci_method,
                    sfs_mode=sfs_mode,
                    workers=num_workers,
                )
                write_output(format_result(result, fmt), output)

                # Generate asymptotic plot if requested
                if plot_asymptotic:
                    from mkado.io.plotting import create_asymptotic_plot

                    try:
                        create_asymptotic_plot(result, plot_asymptotic)
                        typer.echo(f"Asymptotic plot saved to {plot_asymptotic}", err=True)
                    except ValueError as e:
                        typer.echo(f"Could not generate plot: {e}", err=True)
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Aggregated imputed mode
        if use_imputed and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData
            from mkado.analysis.imputed import imputed_mk_test_multi

            worker_results, warnings, _ = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = _collect_gene_data(
                worker_results, "aggregated imputed"
            )

            if gene_data_list:
                result = imputed_mk_test_multi(
                    gene_data=gene_data_list,
                    cutoff=imputed_cutoff,
                    n_bootstrap=bootstrap,
                )
                write_output(format_result(result, fmt), output)
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Per-gene mode
        worker_results, warnings, had_error = run_parallel_batch(
            tasks, num_workers, "Processing alignments"
        )

        for warning in warnings:
            typer.echo(warning, err=True)

        results = [(r.gene_id, r.result) for r in worker_results]

    else:
        # === SEPARATE FILES MODE ===
        validate_polarize_compatibility(
            use_asymptotic=use_asymptotic,
            use_imputed=use_imputed,
            polarize=bool(polarize_pattern),
            polarize_flag="--polarize-pattern",
        )

        ingroup_files = sorted(input_dir.glob(ingroup_pattern))

        if not ingroup_files:
            raise cli_error(f"No files matching '{ingroup_pattern}' in {input_dir}")

        typer.echo(f"Found {len(ingroup_files)} ingroup files", err=True)
        warn_duplicate_stems(ingroup_files)

        # Listed and sorted once: a directory of thousands of genes must not be
        # scanned again for every ingroup file.
        outgroup_candidates = sorted(input_dir.glob(outgroup_pattern))
        outgroup2_candidates = sorted(input_dir.glob(polarize_pattern)) if polarize_pattern else []

        def find_outgroup_file(ingroup_file: Path) -> Path | None:
            base_name = ingroup_file.stem
            if "_ingroup" in base_name:
                outgroup_name = base_name.replace("_ingroup", "_outgroup") + ".fa"
            elif "_in" in base_name:
                outgroup_name = base_name.replace("_in", "_out") + ".fa"
            else:
                outgroup_name = base_name + "_outgroup.fa"

            outgroup_file_path = input_dir / outgroup_name
            if outgroup_file_path.exists():
                return outgroup_file_path
            return find_partner_file(ingroup_file, outgroup_candidates)

        num_workers = get_worker_count(workers, len(ingroup_files))

        tasks = []
        pre_warnings = []
        for ingroup_file in ingroup_files:
            outgroup_file = find_outgroup_file(ingroup_file)
            if outgroup_file is None:
                pre_warnings.append(f"Warning: No outgroup for {ingroup_file.name}")
                continue

            outgroup2_file: Path | None = None
            if polarize_pattern:
                outgroup2_file = find_partner_file(ingroup_file, outgroup2_candidates)
                if outgroup2_file is None:
                    pre_warnings.append(f"Warning: No second outgroup for {ingroup_file.name}")
                    continue

            tasks.append(
                BatchTask(
                    file_path=ingroup_file,
                    outgroup_file=outgroup_file,
                    outgroup2_file=outgroup2_file,
                    reading_frame=reading_frame,
                    use_asymptotic=use_asymptotic,
                    use_imputed=use_imputed and not aggregate,
                    imputed_cutoff=imputed_cutoff,
                    bins=bins,
                    bootstrap=bootstrap,
                    pool_polymorphisms=pool_polymorphisms,
                    min_freq=min_freq,
                    no_singletons=no_singletons,
                    extract_only=is_aggregate,
                    code_table=code_table_id,
                    sfs_mode=sfs_mode,
                    frequency_cutoffs=frequency_cutoffs,
                )
            )

        for warning in pre_warnings:
            typer.echo(warning, err=True)

        if not tasks:
            raise cli_error("No valid file pairs found")

        # Alpha TG mode
        if alpha_tg:
            from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings, _ = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = _collect_gene_data(worker_results, "alpha-TG")

            if gene_data_list:
                result = alpha_tg_from_gene_data(
                    gene_data=gene_data_list,
                    bootstrap_replicates=bootstrap,
                )
                warn_if_alpha_tg_ci_undefined(result)
                write_output(format_result(result, fmt), output)
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Aggregated asymptotic mode
        if use_asymptotic and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings, _ = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = _collect_gene_data(
                worker_results, "aggregated asymptotic"
            )

            if gene_data_list:
                ci_replicates = resolve_ci_replicates(bootstrap, ci_method)
                result = asymptotic_mk_test_aggregated(
                    gene_data=gene_data_list,
                    num_bins=bins,
                    ci_replicates=ci_replicates,
                    frequency_cutoffs=frequency_cutoffs,
                    ci_method=ci_method,
                    sfs_mode=sfs_mode,
                    workers=num_workers,
                )
                write_output(format_result(result, fmt), output)

                # Generate asymptotic plot if requested
                if plot_asymptotic:
                    from mkado.io.plotting import create_asymptotic_plot

                    try:
                        create_asymptotic_plot(result, plot_asymptotic)
                        typer.echo(f"Asymptotic plot saved to {plot_asymptotic}", err=True)
                    except ValueError as e:
                        typer.echo(f"Could not generate plot: {e}", err=True)
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Aggregated imputed mode
        if use_imputed and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData
            from mkado.analysis.imputed import imputed_mk_test_multi

            worker_results, warnings, _ = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = _collect_gene_data(
                worker_results, "aggregated imputed"
            )

            if gene_data_list:
                result = imputed_mk_test_multi(
                    gene_data=gene_data_list,
                    cutoff=imputed_cutoff,
                    n_bootstrap=bootstrap,
                )
                write_output(format_result(result, fmt), output)
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Per-gene mode
        worker_results, warnings, had_error = run_parallel_batch(
            tasks, num_workers, "Processing files"
        )

        for warning in warnings:
            typer.echo(warning, err=True)

        results = [(r.gene_id, r.result) for r in worker_results]

    if results:
        adjusted_pvalues = compute_adjusted_pvalues(results)
        write_batch_output(results, fmt, adjusted_pvalues, output)

        if volcano:
            save_volcano_plot(results, volcano)

        warn_if_plot_asymptotic_ignored(
            plot_asymptotic=plot_asymptotic, use_asymptotic=use_asymptotic
        )
    else:
        report_no_results(had_error)


@app.command()
def info(
    fasta: Annotated[
        Path,
        typer.Argument(help="Path to FASTA file"),
    ],
    reading_frame: Annotated[
        int,
        typer.Option("--reading-frame", "-r", min=1, max=3, help="Reading frame (1-3)"),
    ] = 1,
) -> None:
    """Display information about a FASTA file."""
    from mkado.core.sequences import SequenceSet

    seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)

    typer.echo(f"File: {fasta.name}")
    typer.echo(f"Sequences: {len(seqs)}")

    if seqs.sequences:
        typer.echo(f"Alignment length: {seqs.alignment_length} bp")
        typer.echo(f"Codons: {seqs.num_codons}")
        typer.echo(f"Reading frame: {reading_frame}")

        poly_sites = seqs.polymorphic_codons()
        typer.echo(f"Polymorphic codons: {len(poly_sites)}")

        typer.echo("\nSequences:")
        for seq in seqs.sequences:
            typer.echo(f"  {seq.name} ({len(seq)} bp)")


@app.command()
def codes() -> None:
    """List available NCBI genetic code tables."""
    from mkado.data.genetic_codes import available_code_tables

    typer.echo("Available genetic code tables (use with --code-table):\n")
    for table_id, name, aliases in available_code_tables():
        alias_str = ", ".join(aliases) if aliases else ""
        typer.echo(f"  {table_id:>2}  {name}")
        if alias_str:
            typer.echo(f"      aliases: {alias_str}")


@app.command()
def vcf(
    vcf_file: Annotated[
        Path,
        typer.Option("--vcf", help="Ingroup VCF file (multi-sample, bgzipped+tabix recommended)"),
    ],
    ref: Annotated[
        Path,
        typer.Option(
            "--ref", help="Reference FASTA file, plain or bgzipped (must be faidx-indexed)"
        ),
    ],
    gff: Annotated[
        Path,
        typer.Option("--gff", help="GFF3 annotation file (plain or gzipped)"),
    ],
    outgroup_vcf: Annotated[
        Path,
        typer.Option(
            "--outgroup-vcf",
            help="Single-sample outgroup VCF (same reference; absent positions count as reference)",
        ),
    ],
    # === Gene selection ===
    gene: Annotated[
        Optional[str],
        typer.Option("--gene", help="Single gene ID to analyze"),
    ] = None,
    gene_list: Annotated[
        Optional[Path],
        typer.Option(
            "--gene-list",
            help="File with gene IDs (one per line) to analyze",
            callback=validate_path_not_flag,
        ),
    ] = None,
    # === Analysis type ===
    use_asymptotic: Annotated[
        bool,
        typer.Option("--asymptotic", "-a", help="Use asymptotic MK test"),
    ] = False,
    aggregate: Annotated[
        bool,
        typer.Option(
            "--aggregate/--per-gene",
            help="Aggregate across genes (default) or per-gene results",
        ),
    ] = True,
    alpha_tg: Annotated[
        bool,
        typer.Option("--alpha-tg", help="Compute weighted alpha_TG"),
    ] = False,
    use_imputed: Annotated[
        bool,
        typer.Option("--imputed", help="Use imputed MK test"),
    ] = False,
    # === Asymptotic options ===
    bins: Annotated[
        int,
        typer.Option("--bins", "-b", help="Frequency bins (asymptotic)"),
    ] = 10,
    bootstrap: Annotated[
        int,
        typer.Option("--bootstrap", help="Bootstrap replicates"),
    ] = 100,
    ci_method: CIMethodOption = "monte-carlo",
    sfs_mode: SfsModeOption = "at",
    freq_cutoffs: Annotated[
        str,
        typer.Option("--freq-cutoffs", help="Frequency range 'low,high' (asymptotic)"),
    ] = "0.1,0.9",
    # === Common options ===
    output_format: Annotated[
        str,
        typer.Option("--format", "-f", help="Output format: pretty, tsv, json"),
    ] = "tsv",
    min_freq: Annotated[
        float,
        typer.Option("--min-freq", help="Min derived allele frequency"),
    ] = 0.0,
    no_singletons: Annotated[
        bool,
        typer.Option("--no-singletons", help="Exclude singletons"),
    ] = False,
    code_table: Annotated[
        str,
        typer.Option(
            "--code-table",
            help="Genetic code: name or NCBI table ID",
        ),
    ] = "standard",
    # === Plot options ===
    plot_asymptotic: Annotated[
        Optional[Path],
        typer.Option(
            "--plot-asymptotic",
            help="Generate alpha(x) plot for asymptotic test (PNG, PDF, or SVG)",
            callback=validate_path_not_flag,
        ),
    ] = None,
    volcano: Annotated[
        Optional[Path],
        typer.Option(
            "--volcano",
            help="Generate volcano plot and save to specified path (PNG, PDF, or SVG)",
            callback=validate_path_not_flag,
        ),
    ] = None,
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", min=0, help="Parallel workers (0=auto, 1=sequential)"),
    ] = 0,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", help="Show warnings from htslib/VCF parsing"),
    ] = False,
    output: OutputOption = None,
) -> None:
    """Run MK test from VCF + reference + GFF3 annotation.

    Extracts polymorphism (from ingroup VCF) and divergence (from outgroup VCF)
    data for each gene defined in the GFF3 annotation, then runs the selected
    MK test variant.

    REQUIREMENTS:

    - Ingroup VCF: multi-sample population VCF (bgzipped+tabix recommended)
    - Reference FASTA: genome assembly the VCF was called against (plain or bgzipped, faidx-indexed)
    - GFF3 annotation: gene models with CDS features (plain or gzipped)
    - Outgroup VCF: single-sample VCF of outgroup (for divergence)

    EXAMPLES:

        mkado vcf --vcf pop.vcf.gz --ref genome.fa --gff genes.gff3 --outgroup-vcf outgroup.vcf.gz
        mkado vcf --vcf pop.vcf.gz --ref genome.fa --gff genes.gff3 --outgroup-vcf out.vcf.gz -a
        mkado vcf --vcf pop.vcf.gz --ref genome.fa --gff genes.gff3 --outgroup-vcf out.vcf.gz --gene BRCA1
    """
    # Configure logging — always show warnings (e.g., volcano plot exclusions),
    # --verbose also enables debug messages (e.g., htslib warnings, per-gene details)
    import logging as _logging
    import sys as _sys

    _logging.basicConfig(level=_logging.WARNING, format="%(message)s", stream=_sys.stderr)
    if verbose:
        _logging.getLogger("mkado.io.vcf").setLevel(_logging.DEBUG)
        _logging.getLogger("mkado.io.plotting").setLevel(_logging.DEBUG)

    fmt = resolve_output_format(output_format)

    validate_ci_method(ci_method)
    validate_sfs_mode(sfs_mode)

    # Validate file existence
    for path, name in [(vcf_file, "--vcf"), (ref, "--ref"), (gff, "--gff")]:
        if not path.exists():
            raise cli_error(f"{name} file not found: {path}")

    if not outgroup_vcf.exists():
        raise cli_error(f"--outgroup-vcf file not found: {outgroup_vcf}")

    validate_option_compatibility(
        use_asymptotic=use_asymptotic,
        use_imputed=use_imputed,
        no_singletons=no_singletons,
        min_freq=min_freq,
        alpha_tg=alpha_tg,
    )

    imputed_cutoff = resolve_imputed_cutoff(use_imputed=use_imputed, min_freq=min_freq)

    # Resolve genetic code
    code_table_id = resolve_code_table_or_exit(code_table)

    frequency_cutoffs = parse_frequency_cutoffs(freq_cutoffs)

    # Parse GFF3
    from mkado.io.gff import parse_gff3

    gene_ids = None
    if gene:
        gene_ids = {gene}
    elif gene_list:
        if not gene_list.exists():
            raise cli_error(f"Gene list file not found: {gene_list}")
        gene_ids = set()
        with open(gene_list) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    gene_ids.add(line)
        typer.echo(f"Loaded {len(gene_ids)} gene IDs from {gene_list}", err=True)

    cds_regions = parse_gff3(gff, gene_ids=gene_ids)

    if not cds_regions:
        raise cli_error("No valid CDS regions found in GFF3")

    typer.echo(f"Found {len(cds_regions)} genes in annotation", err=True)

    # Build tasks
    from mkado.vcf_workers import VcfBatchChunk, VcfBatchTask, process_vcf_chunk, process_vcf_gene

    is_aggregate = resolve_is_aggregate(
        use_asymptotic=use_asymptotic,
        use_imputed=use_imputed,
        alpha_tg=alpha_tg,
        aggregate=aggregate,
    )
    warn_if_alpha_tg_ignores_per_gene(alpha_tg=alpha_tg, aggregate=aggregate)

    tasks = [
        VcfBatchTask(
            gene_id=cds.gene_id,
            transcript_id=cds.transcript_id,
            chrom=cds.chrom,
            exons=cds.exons,
            strand=cds.strand,
            phase=cds.phase,
            vcf_path=vcf_file,
            outgroup_vcf_path=outgroup_vcf,
            ref_fasta_path=ref,
            code_table=code_table_id,
            min_freq=min_freq,
            no_singletons=no_singletons,
            use_asymptotic=use_asymptotic and not aggregate,
            bins=bins,
            bootstrap=bootstrap,
            use_imputed=use_imputed and not aggregate,
            imputed_cutoff=imputed_cutoff,
            extract_only=is_aggregate or (gene is None and not use_asymptotic and not use_imputed),
            ci_method=ci_method,
            sfs_mode=sfs_mode,
            frequency_cutoffs=frequency_cutoffs,
        )
        for cds in cds_regions
    ]

    # Determine workers
    num_workers = get_worker_count(workers, len(tasks))

    # Run parallel processing
    worker_results: list[WorkerResult] = []
    batch_warnings: list[str] = []
    had_error = False

    if num_workers == 1:
        with create_rainbow_progress() as progress:
            task_id = progress.add_task("Processing genes", total=len(tasks))
            for task in tasks:
                wr = process_vcf_gene(task)
                if wr.error:
                    batch_warnings.append(wr.error)
                    had_error = True
                elif wr.warning:
                    batch_warnings.append(f"Warning: {wr.warning}")
                if wr.result is not None:
                    worker_results.append(wr)
                progress.advance(task_id)
    else:
        import math
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Partition tasks into chunks (one per worker) so each worker opens
        # VCF/FASTA handles only once and reuses them for all its genes.
        chunk_size = math.ceil(len(tasks) / num_workers)
        chunks = [
            VcfBatchChunk(tasks=tasks[i : i + chunk_size]) for i in range(0, len(tasks), chunk_size)
        ]

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_vcf_chunk, c): c for c in chunks}
            with create_rainbow_progress() as progress:
                task_id = progress.add_task("Processing genes", total=len(tasks))
                for future in as_completed(futures):
                    chunk = futures[future]
                    try:
                        chunk_results = future.result()
                        for wr in chunk_results:
                            if wr.error:
                                batch_warnings.append(wr.error)
                                had_error = True
                            elif wr.warning:
                                batch_warnings.append(f"Warning: {wr.warning}")
                            if wr.result is not None:
                                worker_results.append(wr)
                    except Exception as e:
                        had_error = True
                        for t in chunk.tasks:
                            batch_warnings.append(f"Error processing {t.gene_id}: {e}")
                    progress.advance(task_id, advance=len(chunk.tasks))

    for w in batch_warnings:
        typer.echo(w, err=True)

    if not worker_results:
        report_no_results(had_error)
        return

    # An aggregate mode asked for with --gene still runs below, on that one gene.
    if gene and len(worker_results) == 1 and not is_aggregate:
        write_output(format_result(worker_results[0].result, fmt), output)
        return

    # Alpha TG mode
    if alpha_tg:
        from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
        from mkado.analysis.asymptotic import PolymorphismData

        gene_data_list: list[PolymorphismData] = _collect_gene_data(worker_results, "alpha-TG")
        if gene_data_list:
            result = alpha_tg_from_gene_data(
                gene_data=gene_data_list,
                bootstrap_replicates=bootstrap,
            )
            warn_if_alpha_tg_ci_undefined(result)
            write_output(format_result(result, fmt), output)
        else:
            typer.echo("No valid gene data extracted", err=True)
        return

    # Aggregated asymptotic mode
    if use_asymptotic and aggregate:
        from mkado.analysis.asymptotic import PolymorphismData

        gene_data_list = _collect_gene_data(worker_results, "aggregated asymptotic")
        if gene_data_list:
            ci_replicates = resolve_ci_replicates(bootstrap, ci_method)
            result = asymptotic_mk_test_aggregated(
                gene_data=gene_data_list,
                num_bins=bins,
                ci_replicates=ci_replicates,
                frequency_cutoffs=frequency_cutoffs,
                ci_method=ci_method,
                sfs_mode=sfs_mode,
                workers=num_workers,
            )
            write_output(format_result(result, fmt), output)
            if plot_asymptotic:
                from mkado.io.plotting import create_asymptotic_plot

                try:
                    create_asymptotic_plot(result, plot_asymptotic)
                    typer.echo(f"Asymptotic plot saved to {plot_asymptotic}", err=True)
                except Exception as e:
                    typer.echo(f"Could not generate asymptotic plot: {e}", err=True)
        else:
            typer.echo("No valid gene data extracted", err=True)
        return

    # Aggregated imputed mode
    if use_imputed and aggregate:
        from mkado.analysis.asymptotic import PolymorphismData
        from mkado.analysis.imputed import imputed_mk_test_multi

        gene_data_list = _collect_gene_data(worker_results, "aggregated imputed")
        if gene_data_list:
            result = imputed_mk_test_multi(
                gene_data=gene_data_list,
                cutoff=imputed_cutoff,
                n_bootstrap=bootstrap,
            )
            write_output(format_result(result, fmt), output)
        else:
            typer.echo("No valid gene data extracted", err=True)
        return

    # Per-gene mode: produce standard MK results from PolymorphismData
    from mkado.analysis.asymptotic import PolymorphismData
    from mkado.analysis.mk_test import mk_test_from_counts

    results_list = []
    for wr in worker_results:
        if isinstance(wr.result, PolymorphismData):
            pn = sum(1 for _, t in wr.result.polymorphisms if t == "N")
            ps = sum(1 for _, t in wr.result.polymorphisms if t == "S")
            mk_result = mk_test_from_counts(
                dn=wr.result.dn,
                ds=wr.result.ds,
                pn=pn,
                ps=ps,
            )
            results_list.append((wr.gene_id, mk_result))
        else:
            results_list.append((wr.gene_id, wr.result))

    adjusted_pvalues = compute_adjusted_pvalues(results_list)
    write_batch_output(results_list, fmt, adjusted_pvalues, output)

    if volcano:
        save_volcano_plot(results_list, volcano)

    warn_if_plot_asymptotic_ignored(plot_asymptotic=plot_asymptotic, use_asymptotic=use_asymptotic)


if __name__ == "__main__":
    app()
