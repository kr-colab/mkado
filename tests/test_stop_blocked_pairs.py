"""Tests for surfacing stop-blocked codon pairs (issue #90).

A codon pair with no stop-free mutational ordering is still dropped from
Dn/Ds/Pn/Ps -- that part of the behavior is unchanged. What's new is a
``stop_blocked_pairs`` count, sourced from ``AlignedPair`` (see
test_alignment.py for the counting logic itself) and threaded through every
per-gene result type, mirroring how ``ln``/``ls``/``omega`` were threaded
through in test_omega.py.
"""

from __future__ import annotations

from mkado.analysis.asymptotic import (
    PolymorphismData,
    asymptotic_mk_test,
    asymptotic_mk_test_aggregated,
    extract_polymorphism_data,
)
from mkado.analysis.imputed import imputed_mk_test, imputed_mk_test_multi
from mkado.analysis.mk_test import mk_test
from mkado.analysis.polarized import polarized_mk_test
from mkado.core.codons import GeneticCode

VERTEBRATE_MITO = GeneticCode(table_id=2)


def _write_fasta(tmp_path, name: str, records: list[tuple[str, str]]):
    path = tmp_path / name
    path.write_text("".join(f">{n}\n{s}\n" for n, s in records))
    return path


class TestMKResultStopBlocked:
    def test_full_pipeline(self, tmp_path) -> None:
        ingroup_fa = _write_fasta(tmp_path, "ingroup.fa", [("i1", "AAA")])
        outgroup_fa = _write_fasta(tmp_path, "outgroup.fa", [("o1", "TGG")])

        result = mk_test(ingroup_fa, outgroup_fa, genetic_code=VERTEBRATE_MITO)
        assert (result.dn, result.ds) == (0, 0)
        assert result.stop_blocked_pairs == 1


class TestExtractPolymorphismDataStopBlocked:
    def test_full_pipeline(self, tmp_path) -> None:
        ingroup_fa = _write_fasta(tmp_path, "ingroup.fa", [("i1", "AAA")])
        outgroup_fa = _write_fasta(tmp_path, "outgroup.fa", [("o1", "TGG")])

        data = extract_polymorphism_data(ingroup_fa, outgroup_fa, genetic_code=VERTEBRATE_MITO)
        assert (data.dn, data.ds) == (0, 0)
        assert data.stop_blocked_pairs == 1


class TestAsymptoticStopBlocked:
    def test_single_gene_pipeline(self, tmp_path) -> None:
        ingroup_fa = _write_fasta(tmp_path, "ingroup.fa", [("i1", "AAA")])
        outgroup_fa = _write_fasta(tmp_path, "outgroup.fa", [("o1", "TGG")])

        result = asymptotic_mk_test(ingroup_fa, outgroup_fa, genetic_code=VERTEBRATE_MITO)
        assert (result.dn, result.ds) == (0, 0)
        assert result.stop_blocked_pairs == 1

    def test_aggregated_sums_across_genes(self) -> None:
        gene1 = PolymorphismData(polymorphisms=[], dn=1, ds=1, stop_blocked_pairs=1)
        gene2 = PolymorphismData(polymorphisms=[], dn=1, ds=1, stop_blocked_pairs=2)
        gene3 = PolymorphismData(polymorphisms=[], dn=1, ds=1, stop_blocked_pairs=0)

        result = asymptotic_mk_test_aggregated([gene1, gene2, gene3], num_bins=5)
        assert result.stop_blocked_pairs == 3


class TestPolarizedStopBlocked:
    def test_full_pipeline(self, tmp_path) -> None:
        # outgroup1 and outgroup2 agree on TGG (ancestral), ingroup's AAA is
        # derived: a blocked pair on the ingroup lineage. polarize_fixed_difference
        # and the unpolarized classify_fixed_difference fallback each attempt
        # (and each fail on) the same AAA<->TGG path, so this counts twice.
        ingroup_fa = _write_fasta(tmp_path, "ingroup.fa", [("i1", "AAA")])
        outgroup1_fa = _write_fasta(tmp_path, "outgroup1.fa", [("o1", "TGG")])
        outgroup2_fa = _write_fasta(tmp_path, "outgroup2.fa", [("o2", "TGG")])

        result = polarized_mk_test(
            ingroup_fa, outgroup1_fa, outgroup2_fa, genetic_code=VERTEBRATE_MITO
        )
        assert (result.dn_ingroup, result.ds_ingroup) == (0, 0)
        assert (result.dn_unpolarized, result.ds_unpolarized) == (0, 0)
        assert result.stop_blocked_pairs == 2


class TestImputedStopBlocked:
    def test_single_gene_copies_the_count(self) -> None:
        gene = PolymorphismData(polymorphisms=[], dn=10, ds=5, stop_blocked_pairs=3)
        result = imputed_mk_test(gene)
        assert result.stop_blocked_pairs == 3

    def test_multi_sums_across_genes(self) -> None:
        gene1 = PolymorphismData(polymorphisms=[], dn=1, ds=1, stop_blocked_pairs=2)
        gene2 = PolymorphismData(polymorphisms=[], dn=1, ds=1, stop_blocked_pairs=5)

        result = imputed_mk_test_multi([gene1, gene2])
        assert result.stop_blocked_pairs == 7
