"""Genetic code tables and codon data.

Contains NCBI genetic code tables and codon path computation utilities.
See: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
"""

from __future__ import annotations

from itertools import permutations

# Standard genetic code (NCBI table 1)
STANDARD_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# NCBI genetic code tables
# Each table is defined by its differences from the standard code.
# Format: {table_id: (name, {codon: amino_acid, ...})}
_CODE_DIFFS: dict[int, tuple[str, dict[str, str]]] = {
    2: ("Vertebrate Mitochondrial", {
        "AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W",
    }),
    3: ("Yeast Mitochondrial", {
        "ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T",
        "TGA": "W",
    }),
    4: ("Mold, Protozoan, Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", {
        "TGA": "W",
    }),
    5: ("Invertebrate Mitochondrial", {
        "AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W",
    }),
    6: ("Ciliate, Dasycladacean and Hexamita Nuclear", {
        "TAA": "Q", "TAG": "Q",
    }),
    9: ("Echinoderm and Flatworm Mitochondrial", {
        "AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W",
    }),
    10: ("Euplotid Nuclear", {
        "TGA": "C",
    }),
    11: ("Bacterial, Archaeal and Plant Plastid", {
        # Identical to standard code for amino acid assignments
    }),
    12: ("Alternative Yeast Nuclear", {
        "CTG": "S",
    }),
    13: ("Ascidian Mitochondrial", {
        "AGA": "G", "AGG": "G", "ATA": "M", "TGA": "W",
    }),
    14: ("Alternative Flatworm Mitochondrial", {
        "AAA": "N", "AGA": "S", "AGG": "S", "TAA": "Y", "TGA": "W",
    }),
    16: ("Chlorophycean Mitochondrial", {
        "TAG": "L",
    }),
    21: ("Trematode Mitochondrial", {
        "AAA": "N", "AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W",
    }),
    22: ("Scenedesmus obliquus Mitochondrial", {
        "TAG": "L", "TCA": "*",
    }),
    23: ("Thraustochytrium Mitochondrial", {
        "TTA": "*",
    }),
    24: ("Rhabdopleuridae Mitochondrial", {
        "AGA": "S", "AGG": "K", "TGA": "W",
    }),
    25: ("Candidate Division SR1 and Gracilibacteria", {
        "TGA": "G",
    }),
    26: ("Pachysolen tannophilus Nuclear", {
        "CTG": "A",
    }),
    27: ("Karyorelictea Nuclear", {
        "TAA": "Q", "TAG": "Q", "TGA": "W",
    }),
    29: ("Mesodinium Nuclear", {
        "TAA": "Y", "TAG": "Y",
    }),
    30: ("Peritrich Nuclear", {
        "TAA": "E", "TAG": "E",
    }),
    31: ("Blastocrithidia Nuclear", {
        "TGA": "W",
    }),
    33: ("Cephalodiscidae Mitochondrial UAA-Tyr", {
        "AGA": "S", "AGG": "K", "TAA": "Y", "TGA": "W",
    }),
}


# Short aliases for common genetic codes.
# Keys are lowercase; lookup is case-insensitive.
_CODE_ALIASES: dict[str, int] = {
    "standard": 1,
    "vertebrate-mito": 2,
    "yeast-mito": 3,
    "mold-mito": 4,
    "protozoan-mito": 4,
    "mycoplasma": 4,
    "invertebrate-mito": 5,
    "ciliate": 6,
    "echinoderm-mito": 9,
    "flatworm-mito": 9,
    "euplotid": 10,
    "bacterial": 11,
    "archaeal": 11,
    "plant-plastid": 11,
    "alt-yeast": 12,
    "ascidian-mito": 13,
    "alt-flatworm-mito": 14,
    "chlorophycean-mito": 16,
    "trematode-mito": 21,
    "scenedesmus-mito": 22,
    "thraustochytrium-mito": 23,
    "rhabdopleuridae-mito": 24,
    "sr1": 25,
    "gracilibacteria": 25,
    "pachysolen": 26,
    "karyorelictea": 27,
    "mesodinium": 29,
    "peritrich": 30,
    "blastocrithidia": 31,
    "cephalodiscidae-mito": 33,
}


def resolve_code_table(name_or_id: str) -> int:
    """Resolve a genetic code name or numeric ID to a table ID.

    Accepts either a numeric NCBI table ID (e.g. "2") or a short name
    (e.g. "vertebrate-mito"). Case-insensitive.

    Raises:
        ValueError: If the name or ID is not recognized.
    """
    # Try numeric first
    try:
        table_id = int(name_or_id)
        if table_id == 1 or table_id in _CODE_DIFFS:
            return table_id
        raise ValueError(
            f"Unknown genetic code table {table_id}. "
            f"Run 'mkado codes' to see available tables."
        )
    except ValueError as e:
        if "Unknown genetic code" in str(e):
            raise

    # Try name alias
    key = name_or_id.lower().strip()
    if key in _CODE_ALIASES:
        return _CODE_ALIASES[key]

    raise ValueError(
        f"Unknown genetic code '{name_or_id}'. "
        f"Run 'mkado codes' to see available tables."
    )


def _build_code_table(table_id: int) -> dict[str, str]:
    """Build a full codon table from the standard code plus diffs for the given table ID."""
    if table_id == 1:
        return dict(STANDARD_CODE)
    if table_id not in _CODE_DIFFS:
        msg = (
            f"Unknown genetic code table {table_id}. "
            f"Available tables: {sorted([1] + list(_CODE_DIFFS.keys()))}"
        )
        raise ValueError(msg)
    _, diffs = _CODE_DIFFS[table_id]
    code = dict(STANDARD_CODE)
    code.update(diffs)
    return code


def get_code_table_name(table_id: int) -> str:
    """Get the name of an NCBI genetic code table."""
    if table_id == 1:
        return "Standard"
    if table_id not in _CODE_DIFFS:
        msg = (
            f"Unknown genetic code table {table_id}. "
            f"Available tables: {sorted([1] + list(_CODE_DIFFS.keys()))}"
        )
        raise ValueError(msg)
    return _CODE_DIFFS[table_id][0]


def available_code_tables() -> list[tuple[int, str, list[str]]]:
    """Return list of (table_id, name, aliases) for all available genetic code tables."""
    # Build reverse map: table_id -> list of aliases
    aliases_by_id: dict[int, list[str]] = {}
    for alias, tid in _CODE_ALIASES.items():
        aliases_by_id.setdefault(tid, []).append(alias)

    tables: list[tuple[int, str, list[str]]] = [
        (1, "Standard", aliases_by_id.get(1, []))
    ]
    for table_id in sorted(_CODE_DIFFS.keys()):
        tables.append((table_id, _CODE_DIFFS[table_id][0], aliases_by_id.get(table_id, [])))
    return tables


# All 64 codons in alphabetical order
CODONS = sorted(STANDARD_CODE.keys())

# Build codon lookup table
CODON_TABLE: dict[str, int] = {codon: i for i, codon in enumerate(CODONS)}

# Nucleotides
NUCLEOTIDES = ["A", "C", "G", "T"]


def _compute_codon_paths(
    code: dict[str, str] | None = None,
) -> dict[tuple[str, str], list[tuple[str, int]]]:
    """Compute shortest paths between all codon pairs.

    Args:
        code: Codon-to-amino-acid mapping. Uses standard code if None.

    Returns:
        Dict mapping (codon1, codon2) to a list of (type, position) tuples,
        where type is 'R' for replacement or 'S' for synonymous, and position is
        the codon position (0, 1, or 2) where the change occurs.
    """
    if code is None:
        code = STANDARD_CODE

    paths: dict[tuple[str, str], list[tuple[str, int]]] = {}

    for c1 in CODONS:
        for c2 in CODONS:
            if c1 == c2:
                paths[(c1, c2)] = []
                continue

            # Find positions that differ
            diffs = [i for i in range(3) if c1[i] != c2[i]]

            if len(diffs) == 1:
                # Single nucleotide change
                pos = diffs[0]
                aa1 = code[c1]
                aa2 = code[c2]
                change_type = "S" if aa1 == aa2 else "R"
                paths[(c1, c2)] = [(change_type, pos)]

            elif len(diffs) == 2:
                # Two changes - find shortest path (minimize replacements)
                best_path: list[tuple[str, int]] = []
                best_replacements = 999

                for first_pos in diffs:
                    second_pos = [p for p in diffs if p != first_pos][0]
                    # Make intermediate codon
                    inter = list(c1)
                    inter[first_pos] = c2[first_pos]
                    inter_codon = "".join(inter)

                    if code.get(inter_codon) == "*":
                        continue  # Skip paths through stop codons

                    # Calculate path through this intermediate
                    aa1 = code[c1]
                    aa_inter = code[inter_codon]
                    aa2 = code[c2]

                    step1_type = "S" if aa1 == aa_inter else "R"
                    step2_type = "S" if aa_inter == aa2 else "R"

                    path = [(step1_type, first_pos), (step2_type, second_pos)]
                    n_replacements = sum(1 for t, _ in path if t == "R")

                    if n_replacements < best_replacements:
                        best_replacements = n_replacements
                        best_path = path

                paths[(c1, c2)] = best_path

            else:
                # Three changes - find shortest path
                best_path = []
                best_replacements = 999

                # Try all 6 orderings
                for order in permutations(diffs):
                    current = list(c1)
                    path = []
                    valid = True

                    for i, pos in enumerate(order):
                        old_codon = "".join(current)
                        current[pos] = c2[pos]
                        new_codon = "".join(current)

                        aa_old = code.get(old_codon)
                        aa_new = code.get(new_codon)

                        if aa_new == "*" and i < 2:  # Don't go through stop codons
                            valid = False
                            break

                        change_type = "S" if aa_old == aa_new else "R"
                        path.append((change_type, pos))

                    if valid:
                        n_replacements = sum(1 for t, _ in path if t == "R")
                        if n_replacements < best_replacements:
                            best_replacements = n_replacements
                            best_path = path

                paths[(c1, c2)] = best_path

    return paths


# Pre-computed codon paths for the standard code
CODON_PATHS = _compute_codon_paths()
