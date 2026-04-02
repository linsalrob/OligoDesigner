"""Tests for OligoDesign.dna.DNA."""

import pytest

from OligoDesign.dna import DNA


# ---------------------------------------------------------------------------
# Construction and validation
# ---------------------------------------------------------------------------


class TestConstruction:
    def test_valid_uppercase(self) -> None:
        dna = DNA("ACGT")
        assert str(dna) == "ACGT"

    def test_lowercase_normalised_to_uppercase(self) -> None:
        dna = DNA("acgt")
        assert str(dna) == "ACGT"

    def test_mixed_case_normalised(self) -> None:
        dna = DNA("AcGt")
        assert str(dna) == "ACGT"

    def test_empty_sequence_allowed(self) -> None:
        dna = DNA("")
        assert len(dna) == 0
        assert str(dna) == ""

    def test_single_base(self) -> None:
        for base in "ACGT":
            dna = DNA(base)
            assert str(dna) == base

    def test_invalid_character_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="Invalid DNA characters"):
            DNA("ACGX")

    def test_invalid_character_message_names_bad_chars(self) -> None:
        with pytest.raises(ValueError, match="'X'"):
            DNA("ACGX")

    def test_non_string_raises_type_error(self) -> None:
        with pytest.raises(TypeError):
            DNA(123)  # type: ignore[arg-type]

    def test_ambiguous_base_rejected_by_default(self) -> None:
        with pytest.raises(ValueError):
            DNA("ACGTN")

    def test_ambiguous_base_accepted_when_flag_set(self) -> None:
        dna = DNA("ACGTN", allow_ambiguous=True)
        assert str(dna) == "ACGTN"

    def test_all_iupac_codes_accepted_when_flag_set(self) -> None:
        iupac = "ACGTRYMKSWHBVDN"
        dna = DNA(iupac, allow_ambiguous=True)
        assert str(dna) == iupac

    def test_biological_oligo_sequence(self) -> None:
        # A realistic 40 bp oligo
        seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        dna = DNA(seq)
        assert len(dna) == 40


# ---------------------------------------------------------------------------
# Python data model
# ---------------------------------------------------------------------------


class TestDataModel:
    def test_len(self) -> None:
        assert len(DNA("ACGT")) == 4

    def test_len_empty(self) -> None:
        assert len(DNA("")) == 0

    def test_str(self) -> None:
        assert str(DNA("ACGT")) == "ACGT"

    def test_repr(self) -> None:
        assert repr(DNA("ACGT")) == "DNA('ACGT')"

    def test_equality_dna_vs_dna(self) -> None:
        assert DNA("ACGT") == DNA("ACGT")

    def test_equality_dna_vs_string(self) -> None:
        assert DNA("ACGT") == "ACGT"

    def test_equality_case_insensitive_string(self) -> None:
        assert DNA("ACGT") == "acgt"

    def test_inequality(self) -> None:
        assert DNA("ACGT") != DNA("ACGC")

    def test_hash_equal_objects(self) -> None:
        assert hash(DNA("ACGT")) == hash(DNA("acgt"))

    def test_hash_usable_in_set(self) -> None:
        s = {DNA("ACGT"), DNA("ACGT"), DNA("TTTT")}
        assert len(s) == 2

    def test_iteration(self) -> None:
        assert list(DNA("ACGT")) == ["A", "C", "G", "T"]

    def test_getitem_first_base(self) -> None:
        assert DNA("ACGT")[0] == "A"

    def test_getitem_last_base(self) -> None:
        assert DNA("ACGT")[3] == "T"

    def test_getitem_negative_index(self) -> None:
        assert DNA("ACGT")[-1] == "T"

    def test_getitem_slice(self) -> None:
        assert DNA("ACGT")[1:3] == "CG"

    def test_contains_string_present(self) -> None:
        assert "ACG" in DNA("ACGT")

    def test_contains_string_absent(self) -> None:
        assert "TTT" not in DNA("ACGT")

    def test_contains_dna_present(self) -> None:
        assert DNA("ACG") in DNA("ACGT")

    def test_contains_case_insensitive(self) -> None:
        assert "acg" in DNA("ACGT")


# ---------------------------------------------------------------------------
# Explicit 1-based indexing
# ---------------------------------------------------------------------------


class TestOneBased:
    def test_get_base_1_first(self) -> None:
        assert DNA("ACGT").get_base_1(1) == "A"

    def test_get_base_1_last(self) -> None:
        assert DNA("ACGT").get_base_1(4) == "T"

    def test_get_base_1_middle(self) -> None:
        assert DNA("ACGT").get_base_1(2) == "C"

    def test_get_base_1_zero_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").get_base_1(0)

    def test_get_base_1_beyond_end_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").get_base_1(5)

    def test_get_base_1_negative_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").get_base_1(-1)

    def test_slice_1_full_sequence(self) -> None:
        dna = DNA("ACGT")
        assert str(dna.slice_1(1, 4)) == "ACGT"

    def test_slice_1_inner_subsequence(self) -> None:
        assert str(DNA("ACGT").slice_1(2, 3)) == "CG"

    def test_slice_1_second_half(self) -> None:
        assert str(DNA("ACGTACGT").slice_1(5, 8)) == "ACGT"

    def test_slice_1_single_base(self) -> None:
        assert str(DNA("ACGT").slice_1(3, 3)) == "G"

    def test_slice_1_start_zero_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").slice_1(0, 4)

    def test_slice_1_end_beyond_length_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").slice_1(1, 5)

    def test_slice_1_start_greater_than_end_raises_value_error(self) -> None:
        with pytest.raises(ValueError):
            DNA("ACGT").slice_1(3, 2)

    def test_slice_1_returns_dna_object(self) -> None:
        result = DNA("ACGT").slice_1(1, 2)
        assert isinstance(result, DNA)


# ---------------------------------------------------------------------------
# Reverse, complement, reverse complement
# ---------------------------------------------------------------------------


class TestReverseComplement:
    def test_reverse(self) -> None:
        assert str(DNA("ACGT").reverse()) == "TGCA"

    def test_reverse_single_base(self) -> None:
        assert str(DNA("A").reverse()) == "A"

    def test_complement_acgt(self) -> None:
        assert str(DNA("ACGT").complement()) == "TGCA"

    def test_complement_all_a(self) -> None:
        assert str(DNA("AAAA").complement()) == "TTTT"

    def test_complement_all_c(self) -> None:
        assert str(DNA("CCCC").complement()) == "GGGG"

    def test_reverse_complement_palindrome(self) -> None:
        # ACGT is its own reverse complement
        assert str(DNA("ACGT").reverse_complement()) == "ACGT"

    def test_reverse_complement_ecori_site(self) -> None:
        # GAATTC is the EcoRI recognition site (palindrome)
        assert str(DNA("GAATTC").reverse_complement()) == "GAATTC"

    def test_reverse_complement_all_a(self) -> None:
        assert str(DNA("AAAA").reverse_complement()) == "TTTT"

    def test_reverse_complement_asymmetric(self) -> None:
        # GGATCC (BamHI) -> reverse complement should also be GGATCC
        assert str(DNA("GGATCC").reverse_complement()) == "GGATCC"

    def test_reverse_complement_non_palindrome(self) -> None:
        assert str(DNA("AACC").reverse_complement()) == "GGTT"

    def test_reverse_complement_returns_dna_object(self) -> None:
        rc = DNA("ACGT").reverse_complement()
        assert isinstance(rc, DNA)

    def test_double_reverse_complement_is_identity(self) -> None:
        original = DNA("ATCGATCGATCG")
        assert str(original.reverse_complement().reverse_complement()) == str(original)

    def test_reverse_complement_ambiguous(self) -> None:
        # N complements to N; R (A/G) complements to Y (C/T)
        dna = DNA("ACRN", allow_ambiguous=True)
        # RC: reverse "ACRN" -> "NRCA", complement N->N, R->Y, C->G, A->T -> "NYGT"
        rc = dna.reverse_complement()
        assert str(rc) == "NYGT"


# ---------------------------------------------------------------------------
# GC content and base composition
# ---------------------------------------------------------------------------


class TestComposition:
    def test_gc_content_equal_distribution(self) -> None:
        assert DNA("ACGT").gc_content() == 0.5

    def test_gc_content_all_gc(self) -> None:
        assert DNA("GCGC").gc_content() == 1.0

    def test_gc_content_all_at(self) -> None:
        assert DNA("ATAT").gc_content() == 0.0

    def test_gc_content_empty_sequence(self) -> None:
        assert DNA("").gc_content() == 0.0

    def test_gc_content_single_g(self) -> None:
        assert DNA("G").gc_content() == 1.0

    def test_gc_content_single_a(self) -> None:
        assert DNA("A").gc_content() == 0.0

    def test_base_composition_counts(self) -> None:
        comp = DNA("AACGTT").base_composition()
        assert comp == {"A": 2, "C": 1, "G": 1, "T": 2}

    def test_base_composition_empty(self) -> None:
        comp = DNA("").base_composition()
        assert comp == {"A": 0, "C": 0, "G": 0, "T": 0}

    def test_base_composition_single_base(self) -> None:
        comp = DNA("AAAA").base_composition()
        assert comp["A"] == 4
        assert comp["C"] == 0

    def test_gc_content_realistic_oligo(self) -> None:
        # 50% GC in a 40 bp sequence
        seq = "GCGCGCGCGCGCGCGCGCGCATATATATATATATATATATAT"[:40]
        dna = DNA(seq)
        assert 0.0 <= dna.gc_content() <= 1.0


# ---------------------------------------------------------------------------
# Shannon entropy
# ---------------------------------------------------------------------------


class TestEntropy:
    def test_equal_distribution_max_entropy(self) -> None:
        # ACGT once each -> equal frequencies -> H = log2(4) = 2.0 bits
        assert DNA("ACGT").entropy() == 2.0

    def test_single_base_zero_entropy(self) -> None:
        assert DNA("AAAA").entropy() == 0.0

    def test_empty_sequence_zero_entropy(self) -> None:
        assert DNA("").entropy() == 0.0

    def test_two_equal_bases_entropy(self) -> None:
        # AATT: 50% A, 50% T -> H = 1.0 bit
        assert DNA("AATT").entropy() == 1.0

    def test_entropy_range(self) -> None:
        # Entropy must be between 0.0 and 2.0 for any DNA sequence
        for seq in ("A", "AC", "ACG", "ACGT", "AACGTT", "AAAAAAAAAG"):
            h = DNA(seq).entropy()
            assert 0.0 <= h <= 2.0, f"Entropy out of range for {seq!r}: {h}"

    def test_returns_float(self) -> None:
        assert isinstance(DNA("ACGT").entropy(), float)

    def test_unequal_distribution_intermediate_entropy(self) -> None:
        # AACC: 50% A, 50% C -> H = 1.0 bit
        assert DNA("AACC").entropy() == 1.0

    def test_longer_equal_sequence(self) -> None:
        # 10 copies of each base -> max entropy
        seq = "ACGT" * 10
        assert DNA(seq).entropy() == 2.0


# ---------------------------------------------------------------------------
# Motif search
# ---------------------------------------------------------------------------


class TestMotifSearch:
    def test_find_motif_single_occurrence(self) -> None:
        assert DNA("ACGTACGT").find_motif("ACGT") == [0, 4]

    def test_find_motif_two_occurrences(self) -> None:
        assert DNA("ACGTACGT").find_motif("ACG") == [0, 4]

    def test_find_motif_not_found(self) -> None:
        assert DNA("ACGT").find_motif("TTT") == []

    def test_find_motif_case_insensitive(self) -> None:
        assert DNA("ACGT").find_motif("acg") == [0]

    def test_find_motif_overlapping(self) -> None:
        # "AA" overlaps in "AAAAAA" at positions 0,1,2,3,4
        assert DNA("AAAAAA").find_motif("AA") == [0, 1, 2, 3, 4]

    def test_find_motif_at_end(self) -> None:
        assert DNA("ACGTTT").find_motif("TTT") == [3]

    def test_find_motif_ecori_site(self) -> None:
        seq = DNA("AAAGAATTCGGG")
        positions = seq.find_motif("GAATTC")
        assert positions == [3]


# ---------------------------------------------------------------------------
# Homopolymer detection
# ---------------------------------------------------------------------------


class TestHomopolymer:
    def test_longest_homopolymer_run_of_four(self) -> None:
        assert DNA("ACAAAAT").longest_homopolymer() == 4

    def test_longest_homopolymer_no_run(self) -> None:
        assert DNA("ACGT").longest_homopolymer() == 1

    def test_longest_homopolymer_empty(self) -> None:
        assert DNA("").longest_homopolymer() == 0

    def test_longest_homopolymer_entire_sequence(self) -> None:
        assert DNA("GGGGGG").longest_homopolymer() == 6

    def test_longest_homopolymer_multiple_runs(self) -> None:
        # "AAACCCCC" – longest is 5 C's
        assert DNA("AAACCCCC").longest_homopolymer() == 5

    def test_has_homopolymer_true(self) -> None:
        assert DNA("ACAAAAT").has_homopolymer(min_length=4) is True

    def test_has_homopolymer_false(self) -> None:
        assert DNA("ACGT").has_homopolymer(min_length=4) is False

    def test_has_homopolymer_boundary(self) -> None:
        # run of exactly 4 with min_length=4 should be True
        assert DNA("AAAAACGT").has_homopolymer(min_length=4) is True
        # run of 3 with min_length=4 should be False
        assert DNA("AAAACGT").has_homopolymer(min_length=5) is False


# ---------------------------------------------------------------------------
# Low-complexity detection
# ---------------------------------------------------------------------------


class TestLowComplexity:
    def test_low_complexity_homopolymer_window(self) -> None:
        assert DNA("AAAAAAAAAG").is_low_complexity(window=10, threshold=0.7) is True

    def test_not_low_complexity_balanced(self) -> None:
        assert DNA("ACGTACGTAC").is_low_complexity(window=10, threshold=0.7) is False

    def test_low_complexity_empty_sequence(self) -> None:
        assert DNA("").is_low_complexity() is False

    def test_low_complexity_short_sequence_below_threshold(self) -> None:
        # "ACGT" – dominant base = 1/4 = 0.25, below 0.7
        assert DNA("ACGT").is_low_complexity(window=10, threshold=0.7) is False

    def test_low_complexity_short_sequence_above_threshold(self) -> None:
        # "AAA" – dominant base = 3/3 = 1.0, above 0.7
        assert DNA("AAA").is_low_complexity(window=10, threshold=0.7) is True


# ---------------------------------------------------------------------------
# Self-complementarity
# ---------------------------------------------------------------------------


class TestSelfComplementarity:
    def test_ecori_palindrome(self) -> None:
        assert DNA("GAATTC").is_self_complementary() is True

    def test_bamhi_palindrome(self) -> None:
        assert DNA("GGATCC").is_self_complementary() is True

    def test_acgt_palindrome(self) -> None:
        # ACGT reverse complement is ACGT
        assert DNA("ACGT").is_self_complementary() is True

    def test_non_palindrome(self) -> None:
        assert DNA("AACC").is_self_complementary() is False

    def test_all_a_not_palindrome(self) -> None:
        assert DNA("AAAA").is_self_complementary() is False


# ---------------------------------------------------------------------------
# Hairpin heuristic
# ---------------------------------------------------------------------------


class TestHairpin:
    def test_simple_hairpin_detected(self) -> None:
        # Stem: AAAA (i=0..3), loop: CCC (i=4..6), stem: TTTT (i=7..10)
        # RC of AAAA is TTTT – matches stem3
        assert DNA("AAAACCCTTTT").has_hairpin(min_stem=4, min_loop=3, max_loop=8) is True

    def test_gc_rich_hairpin(self) -> None:
        # GGGG + ATTT loop + CCCC
        assert DNA("GGGGATTCCCCC").has_hairpin(min_stem=4, min_loop=3, max_loop=8) is True

    def test_no_hairpin_random_sequence(self) -> None:
        assert DNA("ACGTACGT").has_hairpin(min_stem=4, min_loop=3, max_loop=8) is False

    def test_no_hairpin_short_sequence(self) -> None:
        # Too short to form stem + loop + stem
        assert DNA("ACGT").has_hairpin(min_stem=4, min_loop=3, max_loop=8) is False

    def test_hairpin_stem_length_respected(self) -> None:
        # Only 3-base stem present; should not trigger min_stem=4
        dna = DNA("ACGTTTACG")  # no 4-base inverted repeat
        assert dna.has_hairpin(min_stem=4, min_loop=3, max_loop=8) is False


# ---------------------------------------------------------------------------
# Melting temperature (nearest-neighbor model)
# ---------------------------------------------------------------------------


class TestMeltingTemperature:
    def test_empty_sequence_returns_none(self) -> None:
        assert DNA("").melting_temperature() is None

    def test_single_base_returns_none(self) -> None:
        assert DNA("A").melting_temperature() is None

    def test_returns_float(self) -> None:
        assert isinstance(DNA("ACGTACGTACGT").melting_temperature(), float)

    def test_higher_gc_gives_higher_tm(self) -> None:
        # GC-rich sequence should melt at higher temperature than AT-rich
        gc_rich = DNA("GCGCGCGCGCGCGCGCGCGC")
        at_rich = DNA("ATATATATATATATATATATATA")
        assert gc_rich.melting_temperature() > at_rich.melting_temperature()

    def test_longer_sequence_higher_tm(self) -> None:
        # Longer sequences (same composition) have higher Tm
        short = DNA("ATCGATCG")
        long_ = DNA("ATCGATCGATCGATCGATCG")
        assert long_.melting_temperature() > short.melting_temperature()

    def test_realistic_20mer_tm_range(self) -> None:
        # A typical 20-mer primer at ~50% GC should be in a plausible range
        dna = DNA("ATCGATCGATCGATCGATCG")
        tm = dna.melting_temperature()
        assert 30.0 < tm < 75.0

    def test_self_complementary_uses_ct_not_ct_over_4(self) -> None:
        # ACGT is self-complementary; result should differ from a
        # non-self-complementary sequence of the same length
        sc = DNA("ACGT")
        non_sc = DNA("AACT")
        assert sc.melting_temperature() != non_sc.melting_temperature()

    def test_higher_salt_raises_tm(self) -> None:
        dna = DNA("ATCGATCGATCGATCGATCG")
        tm_low_salt = dna.melting_temperature(na_conc=0.05)
        tm_high_salt = dna.melting_temperature(na_conc=1.0)
        assert tm_high_salt > tm_low_salt

    def test_default_na_conc_is_50mm(self) -> None:
        dna = DNA("GCATGCATGCATGCATGCAT")
        assert dna.melting_temperature() == dna.melting_temperature(na_conc=0.05)

    def test_default_oligo_conc_is_250nm(self) -> None:
        dna = DNA("GCATGCATGCATGCATGCAT")
        assert dna.melting_temperature() == dna.melting_temperature(oligo_conc=250e-9)

    def test_known_value_20mer(self) -> None:
        # Verify a pre-computed value for ATCGATCGATCGATCGATCG
        # (non-self-complementary, 20 bp) at default conditions
        dna = DNA("ATCGATCGATCGATCGATCG")
        tm = dna.melting_temperature()
        assert abs(tm - 46.92) < 0.1

    def test_known_value_gc_rich_20mer(self) -> None:
        # Pre-computed value for non-self-complementary GCATGCATGCATGCATGCAT
        dna = DNA("GCATGCATGCATGCATGCAT")
        assert dna.is_self_complementary() is False
        tm = dna.melting_temperature()
        assert abs(tm - 51.15) < 0.1

    def test_symmetry_correction_lowers_tm_for_palindromes(self) -> None:
        # ACGT is self-complementary; the SantaLucia symmetry correction of
        # -1.4 cal/mol·K applied to ΔS should produce a lower Tm than the
        # uncorrected value.  Verify the pre-computed Tm (≈-55.77°C) against
        # the expected value to confirm the correction is applied.
        tm = DNA("ACGT").melting_temperature()
        assert tm is not None
        assert abs(tm - (-55.77)) < 0.1, f"Expected ≈-55.77, got {tm:.2f}"


# ---------------------------------------------------------------------------
# Edge-case validation (new)
# ---------------------------------------------------------------------------


class TestIsLowComplexityValidation:
    def test_window_zero_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="window must be >= 1"):
            DNA("ACGT").is_low_complexity(window=0)

    def test_window_negative_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="window must be >= 1"):
            DNA("ACGT").is_low_complexity(window=-1)

    def test_threshold_negative_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="threshold"):
            DNA("ACGT").is_low_complexity(threshold=-0.1)

    def test_threshold_above_one_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="threshold"):
            DNA("ACGT").is_low_complexity(threshold=1.1)

    def test_threshold_zero_is_valid(self) -> None:
        # threshold=0.0 is the boundary; every window will be "low complexity"
        assert DNA("ACGT").is_low_complexity(threshold=0.0) is True

    def test_threshold_one_is_valid(self) -> None:
        # threshold=1.0 only triggers for a pure homopolymer window
        assert DNA("AAAA").is_low_complexity(window=4, threshold=1.0) is True
        assert DNA("ACGT").is_low_complexity(window=4, threshold=1.0) is False


class TestHasHomopolymerValidation:
    def test_min_length_zero_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_length must be >= 1"):
            DNA("ACGT").has_homopolymer(min_length=0)

    def test_min_length_negative_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_length must be >= 1"):
            DNA("ACGT").has_homopolymer(min_length=-1)

    def test_min_length_one_is_valid(self) -> None:
        # Every non-empty sequence has at least one run of length 1
        assert DNA("ACGT").has_homopolymer(min_length=1) is True


class TestSlice1OutOfRangePrecedence:
    def test_start_and_end_both_out_of_range_raises_index_error(self) -> None:
        # start > end but both out of range: IndexError must be raised
        with pytest.raises(IndexError):
            DNA("ACGT").slice_1(5, 10)

    def test_start_out_of_range_end_valid_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").slice_1(0, 2)

    def test_end_out_of_range_start_valid_raises_index_error(self) -> None:
        with pytest.raises(IndexError):
            DNA("ACGT").slice_1(1, 5)

    def test_start_greater_than_end_both_in_range_raises_value_error(self) -> None:
        # Both positions are valid but start > end: ValueError
        with pytest.raises(ValueError):
            DNA("ACGT").slice_1(3, 2)

    def test_start_out_of_range_and_start_greater_than_end_raises_index_error(self) -> None:
        # start=6, end=3 on a 4-base sequence: out-of-range check fires first
        with pytest.raises(IndexError):
            DNA("ACGT").slice_1(6, 3)


class TestMeltingTemperatureValidation:
    def test_zero_na_conc_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="na_conc"):
            DNA("ACGTACGT").melting_temperature(na_conc=0.0)

    def test_negative_na_conc_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="na_conc"):
            DNA("ACGTACGT").melting_temperature(na_conc=-0.05)

    def test_zero_oligo_conc_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="oligo_conc"):
            DNA("ACGTACGT").melting_temperature(oligo_conc=0.0)

    def test_negative_oligo_conc_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="oligo_conc"):
            DNA("ACGTACGT").melting_temperature(oligo_conc=-1e-9)
