"""Tests for OligoDesigner.structured and OligoDesigner.structured_cli."""

from __future__ import annotations

import json
import os
import random

import pytest

from OligoDesigner.structured import (
    StructuredOligo,
    generate_at_rich_palindrome,
    generate_inverted_repeat,
    generate_palindromic_motif,
)
from OligoDesigner.oligo import write_fasta, write_json, write_tsv
from OligoDesigner.structured_cli import main


# ---------------------------------------------------------------------------
# generate_palindromic_motif
# ---------------------------------------------------------------------------


class TestGeneratePalindromicMotif:
    def test_is_palindrome(self) -> None:
        oligo = generate_palindromic_motif(half_length=6, rng=random.Random(1))
        assert oligo.is_palindrome is True

    def test_correct_length_no_spacer(self) -> None:
        oligo = generate_palindromic_motif(half_length=6, spacer_length=0, rng=random.Random(1))
        assert len(oligo.sequence) == 12

    def test_correct_length_with_spacer(self) -> None:
        oligo = generate_palindromic_motif(half_length=6, spacer_length=4, rng=random.Random(1))
        assert len(oligo.sequence) == 16

    def test_oligo_type(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.oligo_type == "palindromic_motif"

    def test_right_arm_is_rc_of_left(self) -> None:
        oligo = generate_palindromic_motif(half_length=8, rng=random.Random(42))
        from OligoDesigner.dna import DNA
        assert oligo.right_arm == str(DNA(oligo.left_arm).reverse_complement())

    def test_spacer_embedded_in_sequence(self) -> None:
        oligo = generate_palindromic_motif(half_length=6, spacer_length=3, rng=random.Random(5))
        assert oligo.sequence == oligo.left_arm + oligo.spacer + oligo.right_arm

    def test_no_spacer_gives_empty_spacer_field(self) -> None:
        oligo = generate_palindromic_motif(half_length=6, spacer_length=0, rng=random.Random(1))
        assert oligo.spacer == ""

    def test_inner_fields_empty(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.inner_left == ""
        assert oligo.inner_right == ""

    def test_only_acgt_bases(self) -> None:
        oligo = generate_palindromic_motif(half_length=10, rng=random.Random(99))
        assert set(oligo.sequence).issubset(set("ACGT"))

    def test_reproducible_with_seed(self) -> None:
        a = generate_palindromic_motif(rng=random.Random(7))
        b = generate_palindromic_motif(rng=random.Random(7))
        assert a.sequence == b.sequence

    def test_half_length_too_small_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_palindromic_motif(half_length=0)

    def test_negative_spacer_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_palindromic_motif(spacer_length=-1)


# ---------------------------------------------------------------------------
# generate_inverted_repeat
# ---------------------------------------------------------------------------


class TestGenerateInvertedRepeat:
    def test_outer_arms_palindromic(self) -> None:
        oligo = generate_inverted_repeat(rng=random.Random(1))
        assert oligo.is_palindrome is True

    def test_inner_core_palindromic(self) -> None:
        oligo = generate_inverted_repeat(rng=random.Random(1))
        assert oligo.inner_is_palindrome is True

    def test_oligo_type(self) -> None:
        oligo = generate_inverted_repeat(rng=random.Random(1))
        assert oligo.oligo_type == "inverted_repeat"

    def test_correct_total_length(self) -> None:
        # outer_arm*2 + inner_half*2 + inner_spacer
        oligo = generate_inverted_repeat(
            inner_half_length=6, outer_arm_length=8, inner_spacer_length=2,
            rng=random.Random(1),
        )
        assert len(oligo.sequence) == 8 + 6 + 2 + 6 + 8

    def test_sequence_structure(self) -> None:
        oligo = generate_inverted_repeat(rng=random.Random(3))
        expected = (
            oligo.left_arm
            + oligo.inner_left
            + oligo.spacer
            + oligo.inner_right
            + oligo.right_arm
        )
        assert oligo.sequence == expected

    def test_right_arm_is_rc_of_left(self) -> None:
        from OligoDesigner.dna import DNA
        oligo = generate_inverted_repeat(outer_arm_length=8, rng=random.Random(10))
        assert oligo.right_arm == str(DNA(oligo.left_arm).reverse_complement())

    def test_inner_right_is_rc_of_inner_left(self) -> None:
        from OligoDesigner.dna import DNA
        oligo = generate_inverted_repeat(inner_half_length=6, rng=random.Random(10))
        assert oligo.inner_right == str(DNA(oligo.inner_left).reverse_complement())

    def test_only_acgt_bases(self) -> None:
        oligo = generate_inverted_repeat(rng=random.Random(77))
        assert set(oligo.sequence).issubset(set("ACGT"))

    def test_reproducible_with_seed(self) -> None:
        a = generate_inverted_repeat(rng=random.Random(55))
        b = generate_inverted_repeat(rng=random.Random(55))
        assert a.sequence == b.sequence

    def test_inner_half_length_too_small_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_inverted_repeat(inner_half_length=0)

    def test_outer_arm_length_too_small_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_inverted_repeat(outer_arm_length=0)

    def test_negative_inner_spacer_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_inverted_repeat(inner_spacer_length=-1)


# ---------------------------------------------------------------------------
# generate_at_rich_palindrome
# ---------------------------------------------------------------------------


class TestGenerateAtRichPalindrome:
    def test_is_palindrome(self) -> None:
        oligo = generate_at_rich_palindrome(rng=random.Random(1))
        assert oligo.is_palindrome is True

    def test_oligo_type(self) -> None:
        oligo = generate_at_rich_palindrome(rng=random.Random(1))
        assert oligo.oligo_type == "at_rich_palindrome"

    def test_left_arm_is_at_only(self) -> None:
        oligo = generate_at_rich_palindrome(half_length=8, rng=random.Random(42))
        assert set(oligo.left_arm).issubset({"A", "T"})

    def test_right_arm_is_at_only(self) -> None:
        oligo = generate_at_rich_palindrome(half_length=8, rng=random.Random(42))
        assert set(oligo.right_arm).issubset({"A", "T"})

    def test_n_spacer_by_default(self) -> None:
        oligo = generate_at_rich_palindrome(spacer_length=3, rng=random.Random(1))
        assert oligo.spacer == "NNN"

    def test_acgt_spacer_when_requested(self) -> None:
        oligo = generate_at_rich_palindrome(
            spacer_length=3, use_n_spacer=False, rng=random.Random(1)
        )
        assert set(oligo.spacer).issubset(set("ACGT"))

    def test_correct_length(self) -> None:
        oligo = generate_at_rich_palindrome(half_length=6, spacer_length=2, rng=random.Random(1))
        assert len(oligo.sequence) == 14

    def test_inner_fields_empty(self) -> None:
        oligo = generate_at_rich_palindrome(rng=random.Random(1))
        assert oligo.inner_left == ""
        assert oligo.inner_right == ""

    def test_sequence_structure(self) -> None:
        oligo = generate_at_rich_palindrome(rng=random.Random(9))
        assert oligo.sequence == oligo.left_arm + oligo.spacer + oligo.right_arm

    def test_reproducible_with_seed(self) -> None:
        a = generate_at_rich_palindrome(rng=random.Random(3))
        b = generate_at_rich_palindrome(rng=random.Random(3))
        assert a.sequence == b.sequence

    def test_half_length_too_small_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_at_rich_palindrome(half_length=0)

    def test_negative_spacer_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_at_rich_palindrome(spacer_length=-1)


# ---------------------------------------------------------------------------
# StructuredOligo serialisation
# ---------------------------------------------------------------------------


class TestStructuredOligoSerialization:
    def _make_set(self) -> list[StructuredOligo]:
        rng = random.Random(0)
        return [
            generate_palindromic_motif(rng=rng),
            generate_inverted_repeat(rng=rng),
            generate_at_rich_palindrome(rng=rng),
        ]

    def test_to_dict_is_dict(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert isinstance(oligo.to_dict(), dict)

    def test_to_dict_has_required_keys(self) -> None:
        d = generate_palindromic_motif(rng=random.Random(1)).to_dict()
        for key in ("name", "sequence", "length", "oligo_type", "is_palindrome", "gc_content", "entropy", "tm"):
            assert key in d

    def test_to_dict_entropy_value(self) -> None:
        d = generate_palindromic_motif(rng=random.Random(1)).to_dict()
        assert isinstance(d["entropy"], float)
        assert 0.0 <= d["entropy"] <= 2.0

    def test_to_dict_tm_is_float(self) -> None:
        d = generate_palindromic_motif(half_length=10, rng=random.Random(1)).to_dict()
        assert isinstance(d["tm"], float)

    def test_tm_realistic_range_for_20bp(self) -> None:
        oligo = generate_palindromic_motif(half_length=10, rng=random.Random(1))
        assert 10.0 < oligo.tm < 85.0

    def test_at_rich_tm_lower_than_gc_rich(self) -> None:
        at_oligo = generate_at_rich_palindrome(half_length=10, rng=random.Random(1))
        gc_oligo = generate_palindromic_motif(half_length=10, rng=random.Random(42))
        # AT-rich palindromes should generally have lower Tm
        # (at_rich_palindrome arms are AT-only)
        assert at_oligo.tm < gc_oligo.tm

    def test_to_tsv_row_is_list_of_strings(self) -> None:
        row = generate_palindromic_motif(rng=random.Random(1)).to_tsv_row()
        assert isinstance(row, list)
        assert all(isinstance(v, str) for v in row)

    def test_tsv_headers_length_matches_row_length(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert len(StructuredOligo.tsv_headers()) == len(oligo.to_tsv_row())

    def test_name_field_in_tsv_row(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        oligo.name = "test_name"
        assert oligo.to_tsv_row()[0] == "test_name"


# ---------------------------------------------------------------------------
# Shared write_fasta / write_json / write_tsv with StructuredOligo
# ---------------------------------------------------------------------------


class TestSharedWriteWithStructuredOligo:
    """Verify that the shared write helpers work with StructuredOligo objects."""

    def _oligos(self) -> list[StructuredOligo]:
        rng = random.Random(0)
        o1 = generate_palindromic_motif(rng=rng)
        o1.name = "p1"
        o2 = generate_at_rich_palindrome(rng=rng)
        o2.name = "a1"
        return [o1, o2]

    def test_write_fasta_creates_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.fa")
        write_fasta(self._oligos(), path)
        assert os.path.exists(path)

    def test_write_fasta_content(self, tmp_path) -> None:
        path = str(tmp_path / "out.fa")
        oligos = self._oligos()
        write_fasta(oligos, path)
        text = open(path).read()
        for oligo in oligos:
            assert f">{oligo.name}\n{oligo.sequence}\n" in text

    def test_write_json_creates_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        write_json(self._oligos(), path)
        assert os.path.exists(path)

    def test_write_json_is_valid_json(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        write_json(self._oligos(), path)
        data = json.loads(open(path).read())
        assert isinstance(data, list)
        assert len(data) == 2

    def test_write_json_has_structured_fields(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        write_json(self._oligos(), path)
        data = json.loads(open(path).read())
        for item in data:
            assert "oligo_type" in item
            assert "left_arm" in item
            assert "is_palindrome" in item
            assert "entropy" in item

    def test_write_tsv_creates_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._oligos(), path)
        assert os.path.exists(path)

    def test_write_tsv_uses_structured_headers(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._oligos(), path)
        lines = open(path).readlines()
        assert lines[0].strip() == "\t".join(StructuredOligo.tsv_headers())

    def test_write_tsv_row_count(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._oligos(), path)
        lines = open(path).readlines()
        assert len(lines) == 3  # header + 2 rows

    def test_write_tsv_empty_list_creates_empty_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv([], path)
        assert os.path.exists(path)
        assert open(path).read() == ""


# ---------------------------------------------------------------------------
# Structured CLI
# ---------------------------------------------------------------------------


class TestStructuredCLI:
    def test_default_run_exits_zero(self) -> None:
        assert main(["--quiet", "--seed", "1"]) == 0

    def test_type_palindrome(self) -> None:
        assert main(["--type", "palindrome", "--count", "3", "--quiet", "--seed", "1"]) == 0

    def test_type_inverted_repeat(self) -> None:
        assert main(["--type", "inverted_repeat", "--count", "3", "--quiet", "--seed", "1"]) == 0

    def test_type_at_rich(self) -> None:
        assert main(["--type", "at_rich", "--count", "3", "--quiet", "--seed", "1"]) == 0

    def test_type_all_produces_3x_count(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main(["--type", "all", "--count", "4", "--json", path, "--quiet", "--seed", "1"])
        data = json.loads(open(path).read())
        assert len(data) == 12  # 4 × 3 types

    def test_fasta_output(self, tmp_path) -> None:
        fasta = str(tmp_path / "out.fa")
        main(["--count", "3", "--type", "palindrome", "--fasta", fasta, "--quiet", "--seed", "1"])
        lines = open(fasta).readlines()
        assert len([l for l in lines if l.startswith(">")])  == 3

    def test_json_output(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main(["--count", "3", "--type", "palindrome", "--json", path, "--quiet", "--seed", "1"])
        data = json.loads(open(path).read())
        assert len(data) == 3
        assert all("oligo_type" in item for item in data)

    def test_tsv_output(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        main(["--count", "3", "--type", "palindrome", "--tsv", path, "--quiet", "--seed", "1"])
        lines = open(path).readlines()
        assert len(lines) == 4  # header + 3 rows

    def test_tsv_uses_structured_headers(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        main(["--count", "2", "--type", "palindrome", "--tsv", path, "--quiet", "--seed", "1"])
        header = open(path).readlines()[0].strip()
        assert header == "\t".join(StructuredOligo.tsv_headers())

    def test_all_outputs_together(self, tmp_path) -> None:
        fasta = str(tmp_path / "out.fa")
        jsn = str(tmp_path / "out.json")
        tsv = str(tmp_path / "out.tsv")
        main([
            "--count", "2", "--type", "all",
            "--fasta", fasta, "--json", jsn, "--tsv", tsv,
            "--quiet", "--seed", "42",
        ])
        assert os.path.exists(fasta)
        assert os.path.exists(jsn)
        assert os.path.exists(tsv)

    def test_reproducible_with_seed(self, tmp_path) -> None:
        f1 = str(tmp_path / "a.fa")
        f2 = str(tmp_path / "b.fa")
        main(["--count", "3", "--fasta", f1, "--quiet", "--seed", "7"])
        main(["--count", "3", "--fasta", f2, "--quiet", "--seed", "7"])
        assert open(f1).read() == open(f2).read()

    def test_invalid_count_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--count", "0"])
        assert exc.value.code != 0

    def test_invalid_half_length_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--half-length", "0"])
        assert exc.value.code != 0

    def test_stdout_summary_printed(self, capsys) -> None:
        main(["--count", "2", "--type", "palindrome", "--seed", "1"])
        captured = capsys.readouterr()
        assert "Generated 2 structured oligos" in captured.out

    def test_prefix_used_in_names(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--count", "2", "--type", "palindrome",
            "--prefix", "myoligo",
            "--json", path,
            "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        assert all(item["name"].startswith("myoligo") for item in data)


# ---------------------------------------------------------------------------
# Edge-case validation: spacer-length enforcement (new)
# ---------------------------------------------------------------------------


class TestSpacerLengthEnforcement:
    """Generators must reject spacer lengths not in SPACER_LENGTHS = (0,2,3,4,5,6)."""

    from OligoDesigner.structured import SPACER_LENGTHS

    def test_palindromic_motif_spacer_1_raises(self) -> None:
        with pytest.raises(ValueError, match="spacer_length"):
            generate_palindromic_motif(spacer_length=1)

    def test_palindromic_motif_spacer_7_raises(self) -> None:
        with pytest.raises(ValueError, match="spacer_length"):
            generate_palindromic_motif(spacer_length=7)

    def test_palindromic_motif_spacer_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="spacer_length"):
            generate_palindromic_motif(spacer_length=-1)

    def test_inverted_repeat_inner_spacer_1_raises(self) -> None:
        with pytest.raises(ValueError, match="inner_spacer_length"):
            generate_inverted_repeat(inner_spacer_length=1)

    def test_inverted_repeat_inner_spacer_7_raises(self) -> None:
        with pytest.raises(ValueError, match="inner_spacer_length"):
            generate_inverted_repeat(inner_spacer_length=7)

    def test_inverted_repeat_inner_spacer_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="inner_spacer_length"):
            generate_inverted_repeat(inner_spacer_length=-1)

    def test_at_rich_palindrome_spacer_1_raises(self) -> None:
        with pytest.raises(ValueError, match="spacer_length"):
            generate_at_rich_palindrome(spacer_length=1)

    def test_at_rich_palindrome_spacer_7_raises(self) -> None:
        with pytest.raises(ValueError, match="spacer_length"):
            generate_at_rich_palindrome(spacer_length=7)

    def test_valid_spacer_0_accepted(self) -> None:
        oligo = generate_palindromic_motif(spacer_length=0, rng=random.Random(1))
        assert oligo.spacer == ""

    def test_valid_spacer_6_accepted(self) -> None:
        oligo = generate_palindromic_motif(spacer_length=6, rng=random.Random(1))
        assert len(oligo.spacer) == 6


class TestCLISpacerLengthEnforcement:
    def test_spacer_length_1_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--spacer-length", "1"])
        assert exc.value.code != 0

    def test_spacer_length_7_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--spacer-length", "7"])
        assert exc.value.code != 0

    def test_spacer_length_valid_0_succeeds(self) -> None:
        assert main(["--spacer-length", "0", "--quiet", "--seed", "1"]) == 0

    def test_spacer_length_valid_4_succeeds(self) -> None:
        assert main(["--spacer-length", "4", "--quiet", "--seed", "1"]) == 0


# ---------------------------------------------------------------------------
# Structured hairpin with N-spacer (new)
# ---------------------------------------------------------------------------


class TestStructuredHairpinNSpacer:
    """Verify that has_hairpin handles N-spacer bases correctly."""

    def test_n_in_loop_does_not_prevent_hairpin_detection(self) -> None:
        # Build an oligo with a clear stem and an N-containing loop.
        # Stem: AAAACCC + N*4 loop + GGGTTTT is NOT a hairpin (no complementarity).
        # Build a real one: AAAACCCC + NNNN + GGGGTTTT should not match.
        # Use a sequence where N occupies loop positions but stems are real.
        # stem5=AAAA loop=NNNNN stem3=TTTT -> has_hairpin must find it
        # since 'N' in stem3 positions means RC('AAAA')='TTTT' != 'NNNN'.
        # So: AAAA + NNN + TTTT should be detected.
        oligo = StructuredOligo(
            sequence="AAAANNNTTTT",
            oligo_type="at_rich_palindrome",
            left_arm="AAAA",
            right_arm="TTTT",
            spacer="NNN",
            inner_left="",
            inner_right="",
        )
        # stem5=AAAA, loop=NNN (3 bases), stem3=TTTT; RC('AAAA')='TTTT' == stem3
        assert oligo.has_hairpin is True

    def test_n_in_stem_position_blocks_hairpin(self) -> None:
        # Replace one stem base with N to break complementarity.
        # Sequence: NAAA + NNN + TTTA  (3-base stem at most, below min_stem=4)
        oligo = StructuredOligo(
            sequence="NAAANNNTTTN",
            oligo_type="at_rich_palindrome",
            left_arm="NAAA",
            right_arm="TTTN",
            spacer="NNN",
            inner_left="",
            inner_right="",
        )
        # No 4-base stem can form because N ≠ complement of any base
        assert oligo.has_hairpin is False

    def test_n_only_spacer_preserves_loop_geometry(self) -> None:
        # generate_at_rich_palindrome with a 4-base spacer should have correct
        # loop geometry.  The N spacer is the loop; the arms are the stem.
        # With half_length=6 the arms are 6 AT bases; the full stem is 6 bp.
        oligo = generate_at_rich_palindrome(
            half_length=6, spacer_length=4, rng=random.Random(0)
        )
        # 6-base AT arm + 4N spacer + 6-base AT arm -> stem len 6, loop len 4
        # has_hairpin should be True (min_stem=4, loop=4 which is in [3,8])
        assert oligo.has_hairpin is True


# ---------------------------------------------------------------------------
# Analysis options: configurable parameters
# ---------------------------------------------------------------------------


class TestAnalysisParameterFields:
    """StructuredOligo analysis parameters are configurable via instance fields."""

    def _make_hairpin_oligo(self) -> StructuredOligo:
        """Return a known-hairpin oligo: AAAA + NNN + TTTT (stem=4, loop=3)."""
        return StructuredOligo(
            sequence="AAAANNNTTTT",
            oligo_type="at_rich_palindrome",
            left_arm="AAAA",
            right_arm="TTTT",
            spacer="NNN",
            inner_left="",
            inner_right="",
        )

    def test_default_min_stem_is_4(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.min_stem == 4

    def test_default_min_loop_is_3(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.min_loop == 3

    def test_default_max_loop_is_8(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.max_loop == 8

    def test_default_min_hp_run_is_4(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.min_hp_run == 4

    def test_high_min_stem_disables_hairpin(self) -> None:
        oligo = self._make_hairpin_oligo()
        assert oligo.has_hairpin is True  # default min_stem=4
        oligo.min_stem = 5  # stem is only 4 bases; raising to 5 disables it
        assert oligo.has_hairpin is False

    def test_high_max_loop_enables_hairpin(self) -> None:
        # Build an oligo with a large loop that is NOT detected with default max_loop=8
        # stem5=AAAA, loop=NNNNNNNNN (9 bases), stem3=TTTT
        oligo = StructuredOligo(
            sequence="AAAANNNNNNNNNTTTT",
            oligo_type="at_rich_palindrome",
            left_arm="AAAA",
            right_arm="TTTT",
            spacer="NNNNNNNNN",
            inner_left="",
            inner_right="",
        )
        assert oligo.has_hairpin is False  # loop=9, above default max_loop=8
        oligo.max_loop = 9
        assert oligo.has_hairpin is True

    def test_has_homopolymer_default(self) -> None:
        # AAAANNNTTTT has a 4-base A run -> has_homopolymer with min_hp_run=4
        oligo = self._make_hairpin_oligo()
        assert oligo.has_homopolymer is True

    def test_has_homopolymer_higher_threshold(self) -> None:
        # Raise min_hp_run to 5: the 4-base A run should no longer be flagged
        oligo = self._make_hairpin_oligo()
        oligo.min_hp_run = 5
        assert oligo.has_homopolymer is False

    def test_complementary_to_default_empty(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        assert oligo.complementary_to == []

    def test_complementary_to_in_to_dict(self) -> None:
        oligo = generate_palindromic_motif(rng=random.Random(1))
        oligo.complementary_to = ["other1"]
        d = oligo.to_dict()
        assert d["complementary_to"] == ["other1"]

    def test_has_homopolymer_in_to_dict(self) -> None:
        oligo = self._make_hairpin_oligo()
        d = oligo.to_dict()
        assert "has_homopolymer" in d
        assert d["has_homopolymer"] is True

    def test_has_homopolymer_in_tsv_headers(self) -> None:
        assert "has_homopolymer" in StructuredOligo.tsv_headers()

    def test_complementary_to_in_tsv_headers(self) -> None:
        assert "complementary_to" in StructuredOligo.tsv_headers()


class TestStructuredCLIAnalysisOptions:
    """Verify that the analysis options are accepted and affect CLI behaviour."""

    def test_min_stem_option_accepted(self) -> None:
        assert main(["--min-stem", "3", "--quiet", "--seed", "1"]) == 0

    def test_min_loop_option_accepted(self) -> None:
        assert main(["--min-loop", "4", "--quiet", "--seed", "1"]) == 0

    def test_max_loop_option_accepted(self) -> None:
        assert main(["--max-loop", "10", "--quiet", "--seed", "1"]) == 0

    def test_min_hp_run_option_accepted(self) -> None:
        assert main(["--min-hp-run", "3", "--quiet", "--seed", "1"]) == 0

    def test_min_overlap_option_accepted(self) -> None:
        assert main(["--min-overlap", "8", "--quiet", "--seed", "1"]) == 0

    def test_high_min_stem_reduces_hairpins(self, tmp_path) -> None:
        path_default = str(tmp_path / "default.json")
        path_high = str(tmp_path / "high.json")

        main(["--count", "10", "--type", "palindrome", "--seed", "7",
              "--json", path_default, "--quiet"])
        main(["--count", "10", "--type", "palindrome", "--seed", "7",
              "--min-stem", "20", "--json", path_high, "--quiet"])

        default_data = json.loads(open(path_default).read())
        high_data = json.loads(open(path_high).read())

        default_hairpins = sum(1 for r in default_data if r["has_hairpin"])
        high_hairpins = sum(1 for r in high_data if r["has_hairpin"])
        # Requiring a 20-base stem on short oligos should yield fewer hairpins
        assert high_hairpins <= default_hairpins

    def test_json_includes_has_homopolymer(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main(["--count", "2", "--type", "palindrome", "--seed", "1",
              "--json", path, "--quiet"])
        data = json.loads(open(path).read())
        for item in data:
            assert "has_homopolymer" in item

    def test_json_includes_complementary_to(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main(["--count", "2", "--type", "palindrome", "--seed", "1",
              "--json", path, "--quiet"])
        data = json.loads(open(path).read())
        for item in data:
            assert "complementary_to" in item

    def test_tsv_includes_has_homopolymer_column(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        main(["--count", "2", "--type", "palindrome", "--seed", "1",
              "--tsv", path, "--quiet"])
        header = open(path).readlines()[0].strip().split("\t")
        assert "has_homopolymer" in header

    def test_tsv_includes_complementary_to_column(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        main(["--count", "2", "--type", "palindrome", "--seed", "1",
              "--tsv", path, "--quiet"])
        header = open(path).readlines()[0].strip().split("\t")
        assert "complementary_to" in header


# ---------------------------------------------------------------------------
# Structured CLI flank / spacer options
# ---------------------------------------------------------------------------


class TestStructuredCLIFlanks:
    """Tests for --five-prime-spacer, --three-prime-spacer, and random-length
    flanking-sequence options on the generate-structured-oligos CLI."""

    def test_five_prime_spacer_prepended(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "3",
            "--five-prime-spacer", "AAAA",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            assert item["sequence"].startswith("AAAA")

    def test_three_prime_spacer_appended(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "3",
            "--three-prime-spacer", "CCCC",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            assert item["sequence"].endswith("CCCC")

    def test_both_spacers_applied(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "2", "--half-length", "5",
            "--spacer-length", "0",
            "--five-prime-spacer", "GGG",
            "--three-prime-spacer", "TTT",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            assert item["sequence"].startswith("GGG")
            assert item["sequence"].endswith("TTT")
            # 2*5 core (no inner spacer) + 3 five-prime + 3 three-prime = 16
            assert item["length"] == 16

    def test_five_prime_random_length(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "3", "--half-length", "5",
            "--spacer-length", "0",
            "--five-prime-random-length", "4",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            # 2*5 core + 4 five-prime = 14
            assert item["length"] == 14

    def test_three_prime_random_length(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "3", "--half-length", "5",
            "--spacer-length", "0",
            "--three-prime-random-length", "6",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            # 2*5 core + 6 three-prime = 16
            assert item["length"] == 16

    def test_random_flanks_are_valid_bases(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "5",
            "--five-prime-random-length", "4",
            "--three-prime-random-length", "4",
            "--json", path, "--quiet", "--seed", "42",
        ])
        data = json.loads(open(path).read())
        valid = set("ACGTN")  # N is allowed inside at_rich palindromes
        for item in data:
            assert set(item["sequence"]).issubset(valid)

    def test_random_flanks_reproducible_with_seed(self, tmp_path) -> None:
        f1 = str(tmp_path / "a.json")
        f2 = str(tmp_path / "b.json")
        main(["--type", "palindrome", "--count", "3", "--five-prime-random-length", "5",
              "--json", f1, "--quiet", "--seed", "99"])
        main(["--type", "palindrome", "--count", "3", "--five-prime-random-length", "5",
              "--json", f2, "--quiet", "--seed", "99"])
        assert open(f1).read() == open(f2).read()

    def test_spacer_lowercase_accepted(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "2",
            "--five-prime-spacer", "acgt",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            assert item["sequence"].startswith("ACGT")

    def test_no_flanks_unchanged(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "2", "--half-length", "5",
            "--spacer-length", "0",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        for item in data:
            assert item["length"] == 10

    def test_mutually_exclusive_five_prime(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "palindrome", "--count", "2",
                "--five-prime-spacer", "AAAA",
                "--five-prime-random-length", "4",
                "--quiet",
            ])
        assert exc.value.code != 0

    def test_mutually_exclusive_three_prime(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "palindrome", "--count", "2",
                "--three-prime-spacer", "TTTT",
                "--three-prime-random-length", "4",
                "--quiet",
            ])
        assert exc.value.code != 0

    def test_invalid_bases_five_prime(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "palindrome", "--count", "2",
                "--five-prime-spacer", "AAAXTT",
                "--quiet",
            ])
        assert exc.value.code != 0

    def test_invalid_bases_three_prime(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "palindrome", "--count", "2",
                "--three-prime-spacer", "NNNN",
                "--quiet",
            ])
        assert exc.value.code != 0

    def test_random_length_zero_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "palindrome", "--count", "2",
                "--five-prime-random-length", "0",
                "--quiet",
            ])
        assert exc.value.code != 0

    # --same-random-oligo tests

    def test_same_random_oligo_all_flanks_equal(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "5",
            "--five-prime-random-length", "4",
            "--same-random-oligo",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        prefixes = [item["sequence"][:4] for item in data]
        assert len(set(prefixes)) == 1, "all oligos should share the same 5' random flank"

    def test_same_random_oligo_three_prime_all_equal(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "5",
            "--three-prime-random-length", "4",
            "--same-random-oligo",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        suffixes = [item["sequence"][-4:] for item in data]
        assert len(set(suffixes)) == 1, "all oligos should share the same 3' random flank"

    def test_without_same_random_oligo_flanks_differ(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "10",
            "--five-prime-random-length", "6",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        prefixes = [item["sequence"][:6] for item in data]
        assert len(set(prefixes)) > 1, "each oligo should get a unique random 5' flank"

    def test_random_flanks_are_all_unique(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "4",
            "--five-prime-random-length", "1",
            "--json", path, "--quiet", "--seed", "1",
        ])
        prefixes = [item["sequence"][:1] for item in json.loads(open(path).read())]
        assert len(set(prefixes)) == 4

    def test_unique_random_flank_capacity_uses_total_count(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "all", "--count", "2",
                "--five-prime-random-length", "1", "--quiet",
            ])
        assert exc.value.code != 0

    def test_same_random_oligo_implies_deduplicate(self, tmp_path, monkeypatch) -> None:
        # Patch the generator in structured_cli's namespace (it was imported by name)
        # so all cores are identical; the shared flank makes every final sequence a duplicate.
        import OligoDesigner.structured_cli as cli_module
        from OligoDesigner import structured as structured_module
        fixed_oligo = structured_module.generate_palindromic_motif(
            half_length=6, rng=random.Random(1)
        )

        def _fixed_palindrome(**kwargs):
            import copy
            return copy.deepcopy(fixed_oligo)

        monkeypatch.setattr(cli_module, "generate_palindromic_motif", _fixed_palindrome)
        path = str(tmp_path / "out.json")
        main([
            "--type", "palindrome", "--count", "3",
            "--five-prime-random-length", "4",
            "--same-random-oligo",
            "--json", path, "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        assert len(data) == 1

    def test_same_random_oligo_requires_random_length(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main([
                "--type", "palindrome", "--count", "2",
                "--five-prime-spacer", "AAAA",
                "--same-random-oligo",
                "--quiet",
            ])
        assert exc.value.code != 0

    def test_same_random_oligo_reproducible_with_seed(self, tmp_path) -> None:
        f1 = str(tmp_path / "a.json")
        f2 = str(tmp_path / "b.json")
        main(["--type", "palindrome", "--count", "4", "--five-prime-random-length", "5",
              "--same-random-oligo", "--json", f1, "--quiet", "--seed", "7"])
        main(["--type", "palindrome", "--count", "4", "--five-prime-random-length", "5",
              "--same-random-oligo", "--json", f2, "--quiet", "--seed", "7"])
        assert open(f1).read() == open(f2).read()
