"""Static data for genetic codes and codon tables."""

from mkado.data.genetic_codes import (
    CODON_TABLE,
    STANDARD_CODE,
    available_code_tables,
    get_code_table_name,
    resolve_code_table,
)

__all__ = [
    "STANDARD_CODE",
    "CODON_TABLE",
    "available_code_tables",
    "get_code_table_name",
    "resolve_code_table",
]
