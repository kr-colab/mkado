"""The output examples in the docs carry the columns the formatters emit today."""

from __future__ import annotations

import re
from pathlib import Path

from mkado.analysis.mk_test import mk_test_from_counts
from mkado.io.output import OutputFormat, format_batch_results, format_result

DOCS = Path(__file__).parent.parent / "docs"


def _header(text: str) -> list[str]:
    return text.splitlines()[0].split("\t")


def test_batch_workflow_lists_every_standard_batch_column() -> None:
    """The Columns sentence and the example block name the batch TSV header in order."""
    result = mk_test_from_counts(dn=1, ds=1, pn=1, ps=1)
    header = _header(format_batch_results([("g", result)], OutputFormat.TSV, [1.0]))
    page = (DOCS / "batch-workflow.rst").read_text()

    sentence = re.search(r"Columns for the standard MK test: (.*?)\. The", page, re.S).group(1)
    assert re.findall(r"``([^`]+)``", sentence) == header

    example = re.search(r"\n   (gene +Dn +Ds .*?)\n", page).group(1)
    assert example.split() == header


def test_dos_page_shows_every_single_result_field() -> None:
    """The TSV and JSON examples on the DoS page carry the single-result fields in order."""
    result = mk_test_from_counts(dn=6, ds=8, pn=1, ps=8)
    page = (DOCS / "dos.rst").read_text()

    example = re.search(r"\n   (Dn +Ds +Pn .*?)\n", page).group(1)
    assert example.split() == _header(format_result(result, OutputFormat.TSV))
    assert re.findall(r'"(\w+)":', page) == list(result.to_dict())
