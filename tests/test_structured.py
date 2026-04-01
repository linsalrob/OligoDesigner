"""Tests for OligoDesign.structured and OligoDesign.structured_cli."""

from __future__ import annotations

import json
import os
import random

import pytest

from OligoDesign.structured import (
    StructuredOligo,
    generate_at_rich_palindrome,
    generate_inverted_repeat,
    generate_palindromic_motif,
)
from OligoDesign.oligo import write_fasta, write_json, write_tsv
from OligoDesign.structured_cli import main


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
        from OligoDesign.dna import DNA
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
        from OligoDesign.dna import DNA
        oligo = generate_inverted_repeat(outer_arm_length=8, rng=random.Random(10))
        assert oligo.right_arm == str(DNA(oligo.left_arm).reverse_complement())

    def test_inner_right_is_rc_of_inner_left(self) -> None:
        from OligoDesign.dna import DNA
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
        for key in ("name", "sequence", "length", "oligo_type", "is_palindrome", "gc_content"):
            assert key in d

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
