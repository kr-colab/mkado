"""Tests for genetic code and codon utilities."""

import pytest

from mkado.core.codons import GeneticCode


class TestGeneticCode:
    """Tests for GeneticCode class."""

    def test_translate_basic(self) -> None:
        """Test basic codon translation."""
        code = GeneticCode()

        assert code.translate("ATG") == "M"  # Start codon
        assert code.translate("TAA") == "*"  # Stop codon
        assert code.translate("TGG") == "W"  # Tryptophan
        assert code.translate("GCT") == "A"  # Alanine

    def test_translate_lowercase(self) -> None:
        """Test that lowercase codons are handled."""
        code = GeneticCode()

        assert code.translate("atg") == "M"

    def test_translate_unknown(self) -> None:
        """Test translation of unknown codons."""
        code = GeneticCode()

        assert code.translate("NNN") == "X"
        assert code.translate("---") == "X"

    def test_translate_sequence(self) -> None:
        """Test translating a full sequence."""
        code = GeneticCode()

        # ATG = M, GCT = A, TAA = *
        seq = "ATGGCTTAA"
        assert code.translate_sequence(seq) == "MA*"

    def test_translate_sequence_reading_frame(self) -> None:
        """Test translation in different reading frames."""
        code = GeneticCode()

        # Frame 1: ATG GCT -> MA
        # Frame 2: TGG CT -> W (incomplete)
        # Frame 3: GGC T -> G (incomplete)
        seq = "ATGGCT"
        assert code.translate_sequence(seq, reading_frame=1) == "MA"
        assert code.translate_sequence(seq, reading_frame=2) == "W"
        assert code.translate_sequence(seq, reading_frame=3) == "G"

    def test_get_path_single_change(self) -> None:
        """Test path for single nucleotide change."""
        code = GeneticCode()

        # AAA (Lys) -> AAG (Lys) - synonymous at position 2
        path = code.get_path("AAA", "AAG")
        assert len(path) == 1
        assert path[0] == ("S", 2)

        # AAA (Lys) -> GAA (Glu) - replacement at position 0
        path = code.get_path("AAA", "GAA")
        assert len(path) == 1
        assert path[0] == ("R", 0)

    def test_get_path_two_changes(self) -> None:
        """Test path for two nucleotide changes."""
        code = GeneticCode()

        # AAA (Lys) -> GAC (Asp), changing positions 0 and 2. Either order passes
        # through a codon of a third amino acid, so both steps are replacements and
        # the two orderings tie. Which one the table keeps is not part of the rule.
        assert sorted(code.get_path("AAA", "GAC")) == [("R", 0), ("R", 2)]

    def test_get_path_picks_the_fewest_replacements(self) -> None:
        """Orderings of a multi-position change are all the same length.

        What separates them is how many steps are replacements, and the table keeps
        the ordering with the fewest. Averaging over the orderings, as Nei-Gojobori
        prescribes for differences, would report more replacements than this.
        """
        code = GeneticCode()

        # AAA (Lys) -> CGC (Arg). Its six orderings carry 1, 2, 3, 3, 3 and 3
        # replacements, so the average is 2.5 and the chosen path reports 1.
        assert code.get_path("AAA", "CGC") == [("R", 1), ("S", 0), ("S", 2)]

    def test_get_path_is_empty_when_every_ordering_hits_a_stop(self) -> None:
        """A pair with no stop-free ordering has no path, so callers drop it.

        Under the vertebrate mitochondrial code this reaches ordinary sense codons:
        AAA is Lys and TGG is Trp, but AGA, AGG and TAA are all stops there.
        """
        code = GeneticCode(table_id=2)

        assert code.translate("AAA") == "K"
        assert code.translate("TGG") == "W"
        assert code.get_path("AAA", "TGG") == []

    def test_get_path_same_codon(self) -> None:
        """Test path for identical codons."""
        code = GeneticCode()

        path = code.get_path("ATG", "ATG")
        assert path == []

    def test_count_synonymous_sites(self) -> None:
        """Test counting synonymous sites."""
        code = GeneticCode()

        # ATG (Met) has no synonymous sites (only one codon for Met)
        assert code.count_synonymous_sites("ATG") == 0.0

        # Leucine codons have variable degeneracy
        # TTT (Phe) - position 2 can be C and still be Phe
        syn_sites = code.count_synonymous_sites("TTT")
        assert syn_sites > 0

    def test_count_synonymous_sites_ambiguous(self) -> None:
        """Test that ambiguous codons return 0 synonymous sites."""
        code = GeneticCode()

        assert code.count_synonymous_sites("NNN") == 0.0
        assert code.count_synonymous_sites("AT-") == 0.0

    def test_is_synonymous_change(self) -> None:
        """Test identifying synonymous changes."""
        code = GeneticCode()

        # TTT -> TTC (both Phe) - synonymous
        assert code.is_synonymous_change("TTT", "TTC") is True

        # TTT -> CTT (Phe -> Leu) - non-synonymous
        assert code.is_synonymous_change("TTT", "CTT") is False

    def test_is_synonymous_change_multiple_diffs(self) -> None:
        """Test that multiple changes return None."""
        code = GeneticCode()

        # Two differences
        assert code.is_synonymous_change("AAA", "GGG") is None


class TestGeneticCodeTables:
    """Tests for alternate genetic code table support."""

    def test_table_id_standard(self) -> None:
        """Standard code via table_id gives same results as default."""
        default = GeneticCode()
        standard = GeneticCode(table_id=1)
        assert default.translate("ATG") == standard.translate("ATG")
        assert default.translate("TGA") == standard.translate("TGA")

    def test_table_id_vertebrate_mito(self) -> None:
        """Vertebrate mitochondrial code has expected differences."""
        code = GeneticCode(table_id=2)
        assert code.translate("AGA") == "*"  # Stop, not Arg
        assert code.translate("AGG") == "*"  # Stop, not Arg
        assert code.translate("ATA") == "M"  # Met, not Ile
        assert code.translate("TGA") == "W"  # Trp, not Stop

    def test_table_id_invertebrate_mito(self) -> None:
        """Invertebrate mitochondrial code has expected differences."""
        code = GeneticCode(table_id=5)
        assert code.translate("AGA") == "S"  # Ser, not Arg
        assert code.translate("TGA") == "W"  # Trp, not Stop

    def test_paths_differ_between_codes(self) -> None:
        """Paths should differ when amino acid assignments differ."""
        std = GeneticCode()
        mito = GeneticCode(table_id=2)
        # ATG->AGA: in standard, AGA=R so this is R; in mito, AGA=* (stop)
        std_path = std.get_path("ATG", "AGA")
        mito_path = mito.get_path("ATG", "AGA")
        assert std_path != mito_path

    def test_codon_paths_cached(self) -> None:
        """Repeated construction with same table_id shares cached paths."""
        code1 = GeneticCode(table_id=5)
        code2 = GeneticCode(table_id=5)
        assert code1._paths is code2._paths

    def test_invalid_table_id(self) -> None:
        """Invalid table ID raises ValueError."""
        with pytest.raises(ValueError, match="Unknown genetic code"):
            GeneticCode(table_id=99)


class TestResolveCodeTable:
    """Tests for resolve_code_table."""

    def test_numeric_id(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        assert resolve_code_table("1") == 1
        assert resolve_code_table("2") == 2
        assert resolve_code_table("5") == 5

    def test_name_alias(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        assert resolve_code_table("standard") == 1
        assert resolve_code_table("vertebrate-mito") == 2
        assert resolve_code_table("invertebrate-mito") == 5

    def test_case_insensitive(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        assert resolve_code_table("Vertebrate-Mito") == 2
        assert resolve_code_table("STANDARD") == 1

    def test_unknown_name_raises(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        with pytest.raises(ValueError, match="Unknown genetic code"):
            resolve_code_table("not-a-code")

    def test_unknown_id_raises(self) -> None:
        from mkado.data.genetic_codes import resolve_code_table

        with pytest.raises(ValueError, match="Unknown genetic code"):
            resolve_code_table("99")


class TestGeneticCodeMemory:
    """Regression tests for #14: GeneticCode instances must not be pinned by their cache.

    The prior `@lru_cache(maxsize=4096)` on the bound `translate` method held a
    strong reference to ``self`` through the cache key, preventing the instance
    from being garbage-collected until the cached method itself was dropped.
    """

    def test_instance_collected_after_translate_use(self) -> None:
        import gc
        import weakref

        code = GeneticCode()
        for codon in ("ATG", "TTT", "TGA", "GCC"):
            code.translate(codon)
            code.count_synonymous_sites(codon)
        ref = weakref.ref(code)
        del code
        gc.collect()
        assert ref() is None, "GeneticCode instance not collected after del + gc.collect()"

    def test_many_instances_collected(self) -> None:
        import gc
        import weakref

        refs = []
        for _ in range(50):
            c = GeneticCode()
            c.translate("ATG")
            c.count_synonymous_sites("TTT")
            refs.append(weakref.ref(c))
            del c
        gc.collect()
        live = sum(1 for r in refs if r() is not None)
        assert live == 0, f"{live}/{len(refs)} GeneticCode instances still alive"
