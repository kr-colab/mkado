"""Errors for problems in the user's input files."""

from __future__ import annotations


class InputError(ValueError):
    """An input file cannot be used as it is.

    The CLI reports these as one error line, not a traceback. They subclass ValueError, so
    code that catches ValueError still catches them.
    """
