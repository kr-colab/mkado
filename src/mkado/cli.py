"""Command-line interface for mkado."""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated, Optional

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
from mkado.batch_workers import BatchTask, WorkerResult, process_gene
from mkado.io.output import OutputFormat, format_batch_results, format_result
from scipy.stats import false_discovery_control

from mkado.analysis.mk_test import MKResult
from mkado.analysis.polarized import PolarizedMKResult

# Console that writes to stderr (so progress doesn't mix with data output)
stderr_console = Console(stderr=True)


def validate_path_not_flag(value: Path | None) -> Path | None:
    """Validate that a Path argument doesn't look like a flag.

    This catches common mistakes like: --option -a (where -a gets consumed as the path)
    """
    if value is not None and str(value).startswith("-"):
        raise typer.BadParameter(
            f"'{value}' looks like a flag, not a file path. Check the order of your arguments."
        )
    return value


def compute_adjusted_pvalues(
    results: list[tuple[str, MKResult | PolarizedMKResult]],
) -> list[float]:
    """Compute Benjamini-Hochberg adjusted p-values for batch results.

    Args:
        results: List of (name, result) tuples from MK tests

    Returns:
        List of adjusted p-values in the same order as input results
    """
    p_values = []
    for _, result in results:
        if isinstance(result, MKResult):
            p_values.append(result.p_value)
        elif isinstance(result, PolarizedMKResult):
            p_values.append(result.p_value_ingroup)
        else:
            p_values.append(1.0)  # Default for unknown types

    if not p_values:
        return []

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


def run_parallel_batch(
    tasks: list[BatchTask],
    num_workers: int,
    description: str,
) -> tuple[list[WorkerResult], list[str]]:
    """Run batch processing with ProcessPoolExecutor."""
    results: list[WorkerResult] = []
    warnings: list[str] = []

    if num_workers == 1:
        with create_rainbow_progress() as progress:
            task_id = progress.add_task(description, total=len(tasks))
            for task in tasks:
                worker_result = process_gene(task)
                if worker_result.error:
                    warnings.append(worker_result.error)
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
                        elif worker_result.warning:
                            warnings.append(f"Warning: {worker_result.warning}")
                        elif worker_result.result is not None:
                            results.append(worker_result)
                    except Exception as e:
                        task = futures[future]
                        warnings.append(f"Error processing {task.file_path.name}: {e}")
                    progress.advance(task_id)

    return results, warnings


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
            help="Exclude singletons (sets --min-freq to 1/n automatically)",
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
        mkado test alignment.fa -i "dmel" -o "dsim" --polarize-match "dyak"
        mkado test ingroup.fa outgroup.fa
        mkado test ingroup.fa outgroup.fa -a
    """
    from mkado.core.sequences import SequenceSet

    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'.", err=True)
        raise typer.Exit(1)

    # Validate option compatibility
    if use_asymptotic and min_freq > 0.0:
        typer.echo(
            "Error: --min-freq cannot be used with --asymptotic. "
            "The asymptotic test uses --freq-cutoffs for frequency filtering.",
            err=True,
        )
        raise typer.Exit(1)

    if use_asymptotic and no_singletons:
        typer.echo(
            "Error: --no-singletons cannot be used with --asymptotic. "
            "The asymptotic test uses --freq-cutoffs for frequency filtering.",
            err=True,
        )
        raise typer.Exit(1)

    if use_imputed and use_asymptotic:
        typer.echo(
            "Error: --imputed and --asymptotic are mutually exclusive.",
            err=True,
        )
        raise typer.Exit(1)

    if use_imputed and no_singletons:
        typer.echo(
            "Error: --no-singletons cannot be used with --imputed. "
            "The imputed test needs low-frequency variants.",
            err=True,
        )
        raise typer.Exit(1)

    if no_singletons and min_freq > 0.0:
        typer.echo(
            "Error: --no-singletons and --min-freq cannot be used together. "
            "--no-singletons automatically sets the frequency threshold.",
            err=True,
        )
        raise typer.Exit(1)

    # Resolve imputed cutoff from --min-freq (default 0.15)
    imputed_cutoff = min_freq if (use_imputed and min_freq > 0.0) else 0.15

    # Build genetic code
    from mkado.core.codons import GeneticCode
    from mkado.data.genetic_codes import resolve_code_table

    try:
        code_table_id = resolve_code_table(code_table)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    genetic_code = GeneticCode(table_id=code_table_id) if code_table_id != 1 else None

    fmt = OutputFormat(output_format)

    # Determine mode
    combined_mode = ingroup_match is not None or outgroup_match is not None

    if combined_mode:
        # Combined file mode
        if not ingroup_match or not outgroup_match:
            typer.echo(
                "Error: Combined mode requires both -i/--ingroup-match and -o/--outgroup-match",
                err=True,
            )
            raise typer.Exit(1)

        if outgroup_file is not None:
            typer.echo(
                "Error: Don't provide outgroup file when using -i/-o (combined mode)",
                err=True,
            )
            raise typer.Exit(1)

        if polarize_file is not None:
            typer.echo(
                "Error: Use --polarize-match instead of -p in combined mode",
                err=True,
            )
            raise typer.Exit(1)

        # Load and filter sequences
        all_seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)
        ingroup_seqs = all_seqs.filter_by_name(ingroup_match)
        outgroup_seqs = all_seqs.filter_by_name(outgroup_match)

        if len(ingroup_seqs) == 0:
            typer.echo(f"Error: No sequences match ingroup pattern '{ingroup_match}'", err=True)
            raise typer.Exit(1)
        if len(outgroup_seqs) == 0:
            typer.echo(f"Error: No sequences match outgroup pattern '{outgroup_match}'", err=True)
            raise typer.Exit(1)

        typer.echo(
            f"Found {len(ingroup_seqs)} ingroup, {len(outgroup_seqs)} outgroup sequences",
            err=True,
        )

        # Calculate singleton frequency threshold if --no-singletons
        if no_singletons:
            n_samples = len(ingroup_seqs)
            if pool_polymorphisms:
                n_samples += len(outgroup_seqs)
            min_freq = 1.0 / n_samples
            typer.echo(
                f"Excluding singletons (min frequency: {min_freq:.4f}, n={n_samples})",
                err=True,
            )

        # Run appropriate test
        if use_asymptotic:
            if polarize_match:
                typer.echo("Error: Polarized asymptotic test not supported", err=True)
                raise typer.Exit(1)
            result = asymptotic_mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                num_bins=bins,
                bootstrap_replicates=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
                genetic_code=genetic_code,
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
            result = imputed_mk_test(poly_data, cutoff=imputed_cutoff)
        elif polarize_match:
            outgroup2_seqs = all_seqs.filter_by_name(polarize_match)
            if len(outgroup2_seqs) == 0:
                typer.echo(
                    f"Error: No sequences match polarize pattern '{polarize_match}'", err=True
                )
                raise typer.Exit(1)
            typer.echo(f"Polarizing with {len(outgroup2_seqs)} outgroup2 sequences", err=True)
            result = polarized_mk_test(
                ingroup=ingroup_seqs,
                outgroup1=outgroup_seqs,
                outgroup2=outgroup2_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                genetic_code=genetic_code,
            )
        else:
            result = mk_test(
                ingroup=ingroup_seqs,
                outgroup=outgroup_seqs,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                genetic_code=genetic_code,
            )

    else:
        # Separate files mode
        if outgroup_file is None:
            typer.echo(
                "Error: Provide outgroup file, or use -i/-o for combined file mode",
                err=True,
            )
            raise typer.Exit(1)

        if polarize_match is not None:
            typer.echo(
                "Error: Use -p/--polarize instead of --polarize-match in separate files mode",
                err=True,
            )
            raise typer.Exit(1)

        # Calculate singleton frequency threshold if --no-singletons
        if no_singletons:
            ingroup_seqs = SequenceSet.from_fasta(fasta, reading_frame=reading_frame)
            outgroup_seqs = SequenceSet.from_fasta(outgroup_file, reading_frame=reading_frame)
            n_samples = len(ingroup_seqs)
            if pool_polymorphisms:
                n_samples += len(outgroup_seqs)
            min_freq = 1.0 / n_samples
            typer.echo(
                f"Excluding singletons (min frequency: {min_freq:.4f}, n={n_samples})",
                err=True,
            )

        # Run appropriate test
        if use_asymptotic:
            if polarize_file:
                typer.echo("Error: Polarized asymptotic test not supported", err=True)
                raise typer.Exit(1)
            result = asymptotic_mk_test(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                num_bins=bins,
                bootstrap_replicates=bootstrap,
                pool_polymorphisms=pool_polymorphisms,
                genetic_code=genetic_code,
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
            result = imputed_mk_test(poly_data, cutoff=imputed_cutoff)
        elif polarize_file:
            result = polarized_mk_test(
                ingroup=fasta,
                outgroup1=outgroup_file,
                outgroup2=polarize_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                genetic_code=genetic_code,
            )
        else:
            result = mk_test(
                ingroup=fasta,
                outgroup=outgroup_file,
                reading_frame=reading_frame,
                pool_polymorphisms=pool_polymorphisms,
                min_frequency=min_freq,
                genetic_code=genetic_code,
            )

    typer.echo(format_result(result, fmt))

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
    elif plot_asymptotic and not use_asymptotic:
        typer.echo("Warning: --plot-asymptotic requires --asymptotic/-a flag", err=True)


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
            help="Exclude singletons (sets frequency threshold to 1/n per gene)",
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
    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'.", err=True)
        raise typer.Exit(1)

    # Resolve genetic code table
    from mkado.data.genetic_codes import resolve_code_table

    try:
        code_table_id = resolve_code_table(code_table)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    # Validate option compatibility
    if use_asymptotic and min_freq > 0.0:
        typer.echo(
            "Error: --min-freq cannot be used with --asymptotic. "
            "The asymptotic test uses --freq-cutoffs for frequency filtering.",
            err=True,
        )
        raise typer.Exit(1)

    if alpha_tg and use_asymptotic:
        typer.echo(
            "Error: --alpha-tg and --asymptotic are mutually exclusive. "
            "Choose one method for estimating alpha.",
            err=True,
        )
        raise typer.Exit(1)

    if use_asymptotic and no_singletons:
        typer.echo(
            "Error: --no-singletons cannot be used with --asymptotic. "
            "The asymptotic test uses --freq-cutoffs for frequency filtering.",
            err=True,
        )
        raise typer.Exit(1)

    if use_imputed and use_asymptotic:
        typer.echo(
            "Error: --imputed and --asymptotic are mutually exclusive.",
            err=True,
        )
        raise typer.Exit(1)

    if use_imputed and alpha_tg:
        typer.echo(
            "Error: --imputed and --alpha-tg are mutually exclusive.",
            err=True,
        )
        raise typer.Exit(1)

    if use_imputed and no_singletons:
        typer.echo(
            "Error: --no-singletons cannot be used with --imputed. "
            "The imputed test needs low-frequency variants.",
            err=True,
        )
        raise typer.Exit(1)

    if no_singletons and min_freq > 0.0:
        typer.echo(
            "Error: --no-singletons and --min-freq cannot be used together. "
            "--no-singletons automatically sets the frequency threshold.",
            err=True,
        )
        raise typer.Exit(1)

    # Resolve imputed cutoff from --min-freq (default 0.15)
    imputed_cutoff = min_freq if (use_imputed and min_freq > 0.0) else 0.15

    fmt = OutputFormat(output_format)

    # Parse frequency cutoffs
    try:
        cutoff_parts = freq_cutoffs.split(",")
        frequency_cutoffs = (float(cutoff_parts[0]), float(cutoff_parts[1]))
    except (ValueError, IndexError):
        typer.echo(f"Error: Invalid frequency cutoffs '{freq_cutoffs}'. Use 'low,high'", err=True)
        raise typer.Exit(1)

    # Auto-detect mode based on -i flag
    combined_mode = ingroup_match is not None

    if combined_mode:
        # === COMBINED FILE MODE ===
        if not outgroup_match:
            typer.echo("Error: -o/--outgroup-match required with -i/--ingroup-match", err=True)
            raise typer.Exit(1)

        if use_asymptotic and polarize_match:
            typer.echo("Error: --asymptotic and --polarize-match are mutually exclusive", err=True)
            raise typer.Exit(1)

        # Find alignment files
        if file_pattern:
            alignment_files = sorted(input_dir.glob(file_pattern))
        else:
            alignment_files = find_alignment_files(input_dir)

        if not alignment_files:
            pattern_msg = file_pattern or "*.fa, *.fasta, *.fna"
            typer.echo(f"No files found in {input_dir} (tried: {pattern_msg})", err=True)
            raise typer.Exit(1)

        typer.echo(f"Found {len(alignment_files)} alignment files", err=True)

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
                extract_only=(use_asymptotic and aggregate)
                or alpha_tg
                or (use_imputed and aggregate),
                code_table=code_table_id,
            )
            for f in alignment_files
        ]

        # Alpha TG mode
        if alpha_tg:
            from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = alpha_tg_from_gene_data(
                    gene_data=gene_data_list,
                    bootstrap_replicates=bootstrap,
                )
                typer.echo(format_result(result, fmt))
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Aggregated asymptotic mode
        if use_asymptotic and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = asymptotic_mk_test_aggregated(
                    gene_data=gene_data_list,
                    num_bins=bins,
                    ci_replicates=bootstrap * 100,
                    frequency_cutoffs=frequency_cutoffs,
                )
                typer.echo(format_result(result, fmt))

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

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = imputed_mk_test_multi(
                    gene_data=gene_data_list,
                    cutoff=imputed_cutoff,
                )
                typer.echo(format_result(result, fmt))
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Per-gene mode
        worker_results, warnings = run_parallel_batch(tasks, num_workers, "Processing alignments")

        for warning in warnings:
            typer.echo(warning, err=True)

        results = [(r.gene_id, r.result) for r in worker_results]

    else:
        # === SEPARATE FILES MODE ===
        if use_asymptotic and polarize_pattern:
            typer.echo(
                "Error: --asymptotic and --polarize-pattern are mutually exclusive", err=True
            )
            raise typer.Exit(1)

        ingroup_files = sorted(input_dir.glob(ingroup_pattern))

        if not ingroup_files:
            typer.echo(f"No files matching '{ingroup_pattern}' in {input_dir}", err=True)
            raise typer.Exit(1)

        typer.echo(f"Found {len(ingroup_files)} ingroup files", err=True)

        def find_outgroup_file(ingroup_file: Path) -> Path | None:
            base_name = ingroup_file.stem
            if "_ingroup" in base_name:
                outgroup_name = base_name.replace("_ingroup", "_outgroup") + ".fa"
            elif "_in" in base_name:
                outgroup_name = base_name.replace("_in", "_out") + ".fa"
            else:
                outgroup_name = base_name + "_outgroup.fa"

            outgroup_file_path = input_dir / outgroup_name
            if not outgroup_file_path.exists():
                potential_matches = list(input_dir.glob(outgroup_pattern))
                prefix = base_name.split("_")[0]
                matches = [f for f in potential_matches if f.stem.startswith(prefix)]
                if matches:
                    return matches[0]
                return None
            return outgroup_file_path

        def find_outgroup2_file(ingroup_file: Path) -> Path | None:
            if not polarize_pattern:
                return None
            base_name = ingroup_file.stem
            potential_matches = list(input_dir.glob(polarize_pattern))
            prefix = base_name.split("_")[0]
            matches = [f for f in potential_matches if f.stem.startswith(prefix)]
            return matches[0] if matches else None

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
                outgroup2_file = find_outgroup2_file(ingroup_file)
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
                    extract_only=(use_asymptotic and aggregate)
                    or alpha_tg
                    or (use_imputed and aggregate),
                    code_table=code_table_id,
                )
            )

        for warning in pre_warnings:
            typer.echo(warning, err=True)

        if not tasks:
            typer.echo("No valid file pairs found", err=True)
            raise typer.Exit(1)

        # Alpha TG mode
        if alpha_tg:
            from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = alpha_tg_from_gene_data(
                    gene_data=gene_data_list,
                    bootstrap_replicates=bootstrap,
                )
                typer.echo(format_result(result, fmt))
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Aggregated asymptotic mode
        if use_asymptotic and aggregate:
            from mkado.analysis.asymptotic import PolymorphismData

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = asymptotic_mk_test_aggregated(
                    gene_data=gene_data_list,
                    num_bins=bins,
                    ci_replicates=bootstrap * 100,
                    frequency_cutoffs=frequency_cutoffs,
                )
                typer.echo(format_result(result, fmt))

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

            worker_results, warnings = run_parallel_batch(
                tasks, num_workers, "Extracting polymorphism data"
            )

            for warning in warnings:
                typer.echo(warning, err=True)

            gene_data_list: list[PolymorphismData] = [
                r.result for r in worker_results if r.result is not None
            ]

            if gene_data_list:
                result = imputed_mk_test_multi(
                    gene_data=gene_data_list,
                    cutoff=imputed_cutoff,
                )
                typer.echo(format_result(result, fmt))
            else:
                typer.echo("No valid gene data extracted", err=True)
            return

        # Per-gene mode
        worker_results, warnings = run_parallel_batch(tasks, num_workers, "Processing files")

        for warning in warnings:
            typer.echo(warning, err=True)

        results = [(r.gene_id, r.result) for r in worker_results]

    if results:
        adjusted_pvalues = compute_adjusted_pvalues(results)
        typer.echo(format_batch_results(results, fmt, adjusted_pvalues))

        # Generate volcano plot if requested
        if volcano:
            from mkado.io.plotting import create_volcano_plot

            try:
                create_volcano_plot(results, volcano)
                typer.echo(f"Volcano plot saved to {volcano}", err=True)
            except ValueError as e:
                typer.echo(f"Could not generate volcano plot: {e}", err=True)

        # Warn if --plot-asymptotic was used without --asymptotic
        if plot_asymptotic and not use_asymptotic:
            typer.echo(
                "Warning: --plot-asymptotic requires --asymptotic/-a flag (ignored)",
                err=True,
            )
    else:
        typer.echo("No results to display", err=True)


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
        typer.Option("--ref", help="Reference FASTA file, plain or bgzipped (must be indexed)"),
    ],
    gff: Annotated[
        Path,
        typer.Option("--gff", help="GFF3 annotation file (plain or gzipped)"),
    ],
    outgroup_vcf: Annotated[
        Path,
        typer.Option(
            "--outgroup-vcf",
            help="Single-sample outgroup VCF (called against same reference)",
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
    workers: Annotated[
        int,
        typer.Option("--workers", "-w", min=0, help="Parallel workers (0=auto, 1=sequential)"),
    ] = 0,
) -> None:
    """Run MK test from VCF + reference + GFF3 annotation.

    Extracts polymorphism (from ingroup VCF) and divergence (from outgroup VCF)
    data for each gene defined in the GFF3 annotation, then runs the selected
    MK test variant.

    REQUIREMENTS:

    - Ingroup VCF: multi-sample population VCF (bgzipped+tabix recommended)
    - Reference FASTA: genome assembly the VCF was called against (plain or bgzipped, indexed)
    - GFF3 annotation: gene models with CDS features (plain or gzipped)
    - Outgroup VCF: single-sample VCF of outgroup (for divergence)

    Install VCF dependencies: pip install mkado[vcf]

    EXAMPLES:

        mkado vcf --vcf pop.vcf.gz --ref genome.fa --gff genes.gff3 --outgroup-vcf outgroup.vcf.gz
        mkado vcf --vcf pop.vcf.gz --ref genome.fa --gff genes.gff3 --outgroup-vcf out.vcf.gz -a
        mkado vcf --vcf pop.vcf.gz --ref genome.fa --gff genes.gff3 --outgroup-vcf out.vcf.gz --gene BRCA1
    """
    # Check for VCF dependencies
    try:
        from mkado.io.vcf import _check_vcf_deps

        _check_vcf_deps()
    except ImportError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    if output_format not in ("pretty", "tsv", "json"):
        typer.echo(f"Error: Invalid format '{output_format}'.", err=True)
        raise typer.Exit(1)

    # Validate file existence
    for path, name in [(vcf_file, "--vcf"), (ref, "--ref"), (gff, "--gff")]:
        if not path.exists():
            typer.echo(f"Error: {name} file not found: {path}", err=True)
            raise typer.Exit(1)

    if not outgroup_vcf.exists():
        typer.echo(f"Error: --outgroup-vcf file not found: {outgroup_vcf}", err=True)
        raise typer.Exit(1)

    # Validate option compatibility (same as batch command)
    if use_asymptotic and min_freq > 0.0:
        typer.echo("Error: --min-freq cannot be used with --asymptotic.", err=True)
        raise typer.Exit(1)

    if use_asymptotic and no_singletons:
        typer.echo("Error: --no-singletons cannot be used with --asymptotic.", err=True)
        raise typer.Exit(1)

    if use_imputed and use_asymptotic:
        typer.echo("Error: --imputed and --asymptotic are mutually exclusive.", err=True)
        raise typer.Exit(1)

    if alpha_tg and use_asymptotic:
        typer.echo("Error: --alpha-tg and --asymptotic are mutually exclusive.", err=True)
        raise typer.Exit(1)

    if use_imputed and alpha_tg:
        typer.echo("Error: --imputed and --alpha-tg are mutually exclusive.", err=True)
        raise typer.Exit(1)

    if use_imputed and no_singletons:
        typer.echo("Error: --no-singletons cannot be used with --imputed.", err=True)
        raise typer.Exit(1)

    if no_singletons and min_freq > 0.0:
        typer.echo("Error: --no-singletons and --min-freq cannot be used together.", err=True)
        raise typer.Exit(1)

    imputed_cutoff = min_freq if (use_imputed and min_freq > 0.0) else 0.15

    # Resolve genetic code
    from mkado.data.genetic_codes import resolve_code_table

    try:
        code_table_id = resolve_code_table(code_table)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    fmt = OutputFormat(output_format)

    # Parse frequency cutoffs
    try:
        cutoff_parts = freq_cutoffs.split(",")
        frequency_cutoffs = (float(cutoff_parts[0]), float(cutoff_parts[1]))
    except (ValueError, IndexError):
        typer.echo(f"Error: Invalid frequency cutoffs '{freq_cutoffs}'.", err=True)
        raise typer.Exit(1)

    # Parse GFF3
    from mkado.io.gff import parse_gff3

    gene_ids = None
    if gene:
        gene_ids = {gene}
    elif gene_list:
        if not gene_list.exists():
            typer.echo(f"Error: Gene list file not found: {gene_list}", err=True)
            raise typer.Exit(1)
        gene_ids = set()
        with open(gene_list) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    gene_ids.add(line)
        typer.echo(f"Loaded {len(gene_ids)} gene IDs from {gene_list}", err=True)

    cds_regions = parse_gff3(gff, gene_ids=gene_ids)

    if not cds_regions:
        typer.echo("Error: No valid CDS regions found in GFF3", err=True)
        raise typer.Exit(1)

    typer.echo(f"Found {len(cds_regions)} genes in annotation", err=True)

    # Build tasks
    from mkado.vcf_workers import VcfBatchTask, process_vcf_gene

    is_aggregate = (use_asymptotic and aggregate) or alpha_tg or (use_imputed and aggregate)

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
        )
        for cds in cds_regions
    ]

    # Single-gene mode
    if gene and len(tasks) == 1:
        tasks[0].extract_only = False
        if not use_asymptotic and not use_imputed:
            tasks[0].extract_only = False

    # Determine workers
    num_workers = get_worker_count(workers, len(tasks))

    # Run parallel processing
    worker_results: list[WorkerResult] = []
    batch_warnings: list[str] = []

    if num_workers == 1:
        with create_rainbow_progress() as progress:
            task_id = progress.add_task("Processing genes", total=len(tasks))
            for task in tasks:
                wr = process_vcf_gene(task)
                if wr.error:
                    batch_warnings.append(wr.error)
                elif wr.warning:
                    batch_warnings.append(f"Warning: {wr.warning}")
                if wr.result is not None:
                    worker_results.append(wr)
                progress.advance(task_id)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_vcf_gene, t): t for t in tasks}
            with create_rainbow_progress() as progress:
                task_id = progress.add_task("Processing genes", total=len(tasks))
                for future in as_completed(futures):
                    try:
                        wr = future.result()
                        if wr.error:
                            batch_warnings.append(wr.error)
                        elif wr.warning:
                            batch_warnings.append(f"Warning: {wr.warning}")
                        if wr.result is not None:
                            worker_results.append(wr)
                    except Exception as e:
                        t = futures[future]
                        batch_warnings.append(f"Error processing {t.gene_id}: {e}")
                    progress.advance(task_id)

    for w in batch_warnings:
        typer.echo(w, err=True)

    if not worker_results:
        typer.echo("No results to display", err=True)
        raise typer.Exit(1)

    # Single gene mode: output single result
    if gene and len(worker_results) == 1:
        typer.echo(format_result(worker_results[0].result, fmt))
        return

    # Alpha TG mode
    if alpha_tg:
        from mkado.analysis.alpha_tg import alpha_tg_from_gene_data
        from mkado.analysis.asymptotic import PolymorphismData

        gene_data_list: list[PolymorphismData] = [
            r.result for r in worker_results if r.result is not None
        ]
        if gene_data_list:
            result = alpha_tg_from_gene_data(
                gene_data=gene_data_list,
                bootstrap_replicates=bootstrap,
            )
            typer.echo(format_result(result, fmt))
        else:
            typer.echo("No valid gene data extracted", err=True)
        return

    # Aggregated asymptotic mode
    if use_asymptotic and aggregate:
        from mkado.analysis.asymptotic import PolymorphismData

        gene_data_list = [r.result for r in worker_results if r.result is not None]
        if gene_data_list:
            result = asymptotic_mk_test_aggregated(
                gene_data=gene_data_list,
                num_bins=bins,
                ci_replicates=bootstrap * 100,
                frequency_cutoffs=frequency_cutoffs,
            )
            typer.echo(format_result(result, fmt))
        else:
            typer.echo("No valid gene data extracted", err=True)
        return

    # Aggregated imputed mode
    if use_imputed and aggregate:
        from mkado.analysis.asymptotic import PolymorphismData
        from mkado.analysis.imputed import imputed_mk_test_multi

        gene_data_list = [r.result for r in worker_results if r.result is not None]
        if gene_data_list:
            result = imputed_mk_test_multi(
                gene_data=gene_data_list,
                cutoff=imputed_cutoff,
            )
            typer.echo(format_result(result, fmt))
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

    if results_list:
        adjusted_pvalues = compute_adjusted_pvalues(results_list)
        typer.echo(format_batch_results(results_list, fmt, adjusted_pvalues))
    else:
        typer.echo("No results to display", err=True)


if __name__ == "__main__":
    app()
