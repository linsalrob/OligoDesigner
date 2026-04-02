"""Tests for OligoDesign.oligo and OligoDesign.cli."""

from __future__ import annotations

import json
import os
import random

import pytest

from OligoDesign.dna import DNA
from OligoDesign.oligo import (
    OligoAnalysis,
    analyse_oligo,
    find_complementary_pairs,
    has_tandem_repeat,
    random_oligo,
    read_json,
    write_fasta,
    write_json,
    write_tsv,
)
from OligoDesign.cli import main


# ---------------------------------------------------------------------------
# random_oligo
# ---------------------------------------------------------------------------


class TestRandomOligo:
    def test_default_length_is_40(self) -> None:
        oligo = random_oligo(rng=random.Random(1))
        assert len(oligo) == 40

    def test_custom_length(self) -> None:
        oligo = random_oligo(length=20, rng=random.Random(1))
        assert len(oligo) == 20

    def test_zero_length(self) -> None:
        oligo = random_oligo(length=0, rng=random.Random(1))
        assert len(oligo) == 0

    def test_returns_dna_instance(self) -> None:
        oligo = random_oligo(rng=random.Random(1))
        assert isinstance(oligo, DNA)

    def test_only_valid_bases(self) -> None:
        oligo = random_oligo(length=100, rng=random.Random(42))
        assert set(str(oligo)).issubset({"A", "C", "G", "T"})

    def test_reproducible_with_seed(self) -> None:
        a = random_oligo(rng=random.Random(99))
        b = random_oligo(rng=random.Random(99))
        assert str(a) == str(b)

    def test_different_seeds_differ(self) -> None:
        a = random_oligo(length=100, rng=random.Random(1))
        b = random_oligo(length=100, rng=random.Random(2))
        assert str(a) != str(b)

    def test_negative_length_raises(self) -> None:
        with pytest.raises(ValueError):
            random_oligo(length=-1)


# ---------------------------------------------------------------------------
# has_tandem_repeat
# ---------------------------------------------------------------------------


class TestHasTandemRepeat:
    def test_dinucleotide_repeat_detected(self) -> None:
        # ATATAT = AT * 3
        assert has_tandem_repeat(DNA("ATATAT")) is True

    def test_trinucleotide_repeat_detected(self) -> None:
        # ATCATCATC = ATC * 3
        assert has_tandem_repeat(DNA("ATCATCATC")) is True

    def test_no_tandem_repeat(self) -> None:
        # ACGTACGT only 2 copies of ACGT, below min_count=3
        assert has_tandem_repeat(DNA("ACGTACGT")) is False

    def test_short_sequence_no_repeat(self) -> None:
        assert has_tandem_repeat(DNA("ACGT")) is False

    def test_repeat_embedded_in_longer_sequence(self) -> None:
        # GGG + ATCATCATC + TTT
        assert has_tandem_repeat(DNA("GGGATCATCATCTTT")) is True

    def test_min_count_param(self) -> None:
        # ACGTACGT is 2 copies; should be found with min_count=2
        assert has_tandem_repeat(DNA("ACGTACGT"), min_unit=4, max_unit=4, min_count=2) is True

    def test_exact_minimum_copies(self) -> None:
        # ATATAT has exactly 3 copies of AT (the default min_count)
        assert has_tandem_repeat(DNA("ATATAT"), min_unit=2, max_unit=2, min_count=3) is True


# ---------------------------------------------------------------------------
# find_complementary_pairs
# ---------------------------------------------------------------------------


class TestFindComplementaryPairs:
    def test_perfect_complement_detected(self) -> None:
        # oligo1 is the reverse complement of oligo2 (20 bp overlap)
        seq = "ATCGATCGATCGATCGATCG"
        rc = str(DNA(seq).reverse_complement())
        oligos = [DNA(seq), DNA(rc)]
        names = ["o1", "o2"]
        pairs = find_complementary_pairs(oligos, names, min_overlap=10)
        assert "o2" in pairs["o1"]
        assert "o1" in pairs["o2"]

    def test_unrelated_oligos_no_pair(self) -> None:
        # Use short sequences below min_overlap so complementarity window fails
        oligos = [DNA("ACGTACGT"), DNA("TTTTTTTT")]
        names = ["a", "b"]
        pairs = find_complementary_pairs(oligos, names, min_overlap=10)
        assert pairs["a"] == []
        assert pairs["b"] == []

    def test_returns_dict_with_all_names(self) -> None:
        oligos = [DNA("ACGT"), DNA("TTTT"), DNA("GGGG")]
        names = ["x", "y", "z"]
        pairs = find_complementary_pairs(oligos, names, min_overlap=1)
        assert set(pairs.keys()) == {"x", "y", "z"}

    def test_internal_complementarity_detected(self) -> None:
        # The complementary region (ATCG / CGAT) is buried inside both oligos,
        # not at either 5' or 3' end.  The function must detect it.
        seq_a = "AAAAATCGAAAAAAAAA"   # ATCG at positions 5-8
        seq_b = "GGGGGCGATGGGGGGGGG"  # CGAT at positions 5-8 (RC of ATCG)
        oligos = [DNA(seq_a), DNA(seq_b)]
        names = ["o1", "o2"]
        pairs = find_complementary_pairs(oligos, names, min_overlap=4)
        assert "o2" in pairs["o1"]
        assert "o1" in pairs["o2"]

    def test_short_oligos_below_min_overlap_no_pair(self) -> None:
        # Both sequences are shorter than min_overlap; no pair should be reported
        oligos = [DNA("ACGT"), DNA("ACGT")]
        names = ["a", "b"]
        pairs = find_complementary_pairs(oligos, names, min_overlap=10)
        assert pairs["a"] == []
        assert pairs["b"] == []


# ---------------------------------------------------------------------------
# analyse_oligo
# ---------------------------------------------------------------------------


class TestAnalyseOligo:
    def test_returns_oligo_analysis(self) -> None:
        result = analyse_oligo(DNA("ACGT"), name="test")
        assert isinstance(result, OligoAnalysis)

    def test_name_propagated(self) -> None:
        result = analyse_oligo(DNA("ACGT"), name="my_oligo")
        assert result.name == "my_oligo"

    def test_sequence_propagated(self) -> None:
        result = analyse_oligo(DNA("acgt"), name="t")
        assert result.sequence == "ACGT"

    def test_length_correct(self) -> None:
        result = analyse_oligo(DNA("ACGTACGT"), name="t")
        assert result.length == 8

    def test_gc_content_correct(self) -> None:
        result = analyse_oligo(DNA("ACGT"), name="t")
        assert result.gc_content == 0.5

    def test_palindrome_detected(self) -> None:
        result = analyse_oligo(DNA("GAATTC"), name="ecori")
        assert result.is_palindrome is True

    def test_non_palindrome(self) -> None:
        result = analyse_oligo(DNA("AACC"), name="t")
        assert result.is_palindrome is False

    def test_homopolymer_detected(self) -> None:
        result = analyse_oligo(DNA("ACAAAAGT"), name="t")
        assert result.has_homopolymer is True

    def test_hairpin_detected(self) -> None:
        result = analyse_oligo(DNA("AAAACCCTTTT"), name="t")
        assert result.has_hairpin is True

    def test_tandem_repeat_detected(self) -> None:
        result = analyse_oligo(DNA("ATATAT"), name="t")
        assert result.has_tandem_repeat is True

    def test_complementary_to_empty_by_default(self) -> None:
        result = analyse_oligo(DNA("ACGT"), name="t")
        assert result.complementary_to == []

    def test_base_composition_correct(self) -> None:
        result = analyse_oligo(DNA("AACGTT"), name="t")
        assert result.base_composition == {"A": 2, "C": 1, "G": 1, "T": 2}

    def test_entropy_present(self) -> None:
        result = analyse_oligo(DNA("ACGT"), name="t")
        assert isinstance(result.entropy, float)

    def test_entropy_max_for_equal_bases(self) -> None:
        result = analyse_oligo(DNA("ACGT"), name="t")
        assert result.entropy == 2.0

    def test_entropy_zero_for_homopolymer(self) -> None:
        result = analyse_oligo(DNA("AAAA"), name="t")
        assert result.entropy == 0.0

    def test_tm_present(self) -> None:
        result = analyse_oligo(DNA("ATCGATCGATCGATCGATCG"), name="t")
        assert isinstance(result.tm, float)

    def test_tm_realistic_range_for_20mer(self) -> None:
        result = analyse_oligo(DNA("ATCGATCGATCGATCGATCG"), name="t")
        assert 30.0 < result.tm < 75.0


# ---------------------------------------------------------------------------
# OligoAnalysis serialisation
# ---------------------------------------------------------------------------


class TestOligoAnalysisSerialization:
    def _make(self) -> OligoAnalysis:
        return analyse_oligo(DNA("GAATTC"), name="ecori")

    def test_to_dict_is_dict(self) -> None:
        assert isinstance(self._make().to_dict(), dict)

    def test_to_dict_has_expected_keys(self) -> None:
        d = self._make().to_dict()
        for key in ("name", "sequence", "length", "gc_content", "entropy", "tm", "is_palindrome"):
            assert key in d

    def test_to_dict_entropy_value(self) -> None:
        d = self._make().to_dict()
        assert isinstance(d["entropy"], float)
        assert 0.0 <= d["entropy"] <= 2.0

    def test_to_dict_tm_is_float(self) -> None:
        d = self._make().to_dict()
        assert isinstance(d["tm"], float)

    def test_to_tsv_row_is_list_of_strings(self) -> None:
        row = self._make().to_tsv_row()
        assert isinstance(row, list)
        assert all(isinstance(v, str) for v in row)

    def test_tsv_headers_length_matches_row_length(self) -> None:
        a = self._make()
        assert len(OligoAnalysis.tsv_headers()) == len(a.to_tsv_row())


# ---------------------------------------------------------------------------
# write_fasta / write_json / write_tsv
# ---------------------------------------------------------------------------


class TestWriteOutputs:
    def _analyses(self) -> list[OligoAnalysis]:
        return [
            analyse_oligo(DNA("ACGT"), name="o1"),
            analyse_oligo(DNA("TTTT"), name="o2"),
        ]

    def test_write_fasta_creates_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.fa")
        write_fasta(self._analyses(), path)
        assert os.path.exists(path)

    def test_write_fasta_content(self, tmp_path) -> None:
        path = str(tmp_path / "out.fa")
        write_fasta(self._analyses(), path)
        text = open(path).read()
        assert ">o1\nACGT\n" in text
        assert ">o2\nTTTT\n" in text

    def test_write_json_creates_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        write_json(self._analyses(), path)
        assert os.path.exists(path)

    def test_write_json_is_valid_json(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        write_json(self._analyses(), path)
        data = json.loads(open(path).read())
        assert isinstance(data, list)
        assert len(data) == 2

    def test_write_json_has_required_fields(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        write_json(self._analyses(), path)
        data = json.loads(open(path).read())
        for item in data:
            for key in ("name", "sequence", "length", "gc_content", "entropy", "tm"):
                assert key in item

    def test_write_json_tm_is_null_for_short_sequences(self, tmp_path) -> None:
        # Single-base and 4-base sequences are borderline; ensure JSON output
        # is valid RFC 8259 (no bare NaN).
        path = str(tmp_path / "out.json")
        write_json(self._analyses(), path)
        text = open(path).read()
        assert "NaN" not in text
        assert "Infinity" not in text

    def test_write_tsv_creates_file(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._analyses(), path)
        assert os.path.exists(path)

    def test_write_tsv_has_header(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._analyses(), path)
        lines = open(path).readlines()
        assert lines[0].strip() == "\t".join(OligoAnalysis.tsv_headers())

    def test_write_tsv_row_count(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._analyses(), path)
        lines = open(path).readlines()
        # header + 2 data rows
        assert len(lines) == 3

    def test_write_tsv_tab_separated(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        write_tsv(self._analyses(), path)
        lines = open(path).readlines()
        assert "\t" in lines[1]

    def test_write_tsv_mixed_schema_raises(self, tmp_path) -> None:
        """Mixing OligoAnalysis with StructuredOligo in one write_tsv call must fail."""
        from OligoDesign.structured import generate_palindromic_motif
        import random as _random
        path = str(tmp_path / "out.tsv")
        mixed = [
            analyse_oligo(DNA("ACGT"), name="o1"),
            generate_palindromic_motif(rng=_random.Random(1)),
        ]
        with pytest.raises(TypeError, match="Mixed TSV schemas"):
            write_tsv(mixed, path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestCLI:
    def test_default_run_exits_zero(self) -> None:
        assert main(["--quiet", "--seed", "1"]) == 0

    def test_custom_count_and_length(self) -> None:
        assert main(["--count", "5", "--length", "20", "--quiet", "--seed", "1"]) == 0

    def test_fasta_output(self, tmp_path) -> None:
        fasta = str(tmp_path / "out.fa")
        main(["--count", "3", "--length", "20", "--fasta", fasta, "--quiet", "--seed", "1"])
        lines = open(fasta).readlines()
        header_lines = [l for l in lines if l.startswith(">")]
        assert len(header_lines) == 3

    def test_json_output(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main(["--count", "3", "--length", "20", "--json", path, "--quiet", "--seed", "1"])
        data = json.loads(open(path).read())
        assert len(data) == 3

    def test_tsv_output(self, tmp_path) -> None:
        path = str(tmp_path / "out.tsv")
        main(["--count", "3", "--length", "20", "--tsv", path, "--quiet", "--seed", "1"])
        lines = open(path).readlines()
        assert len(lines) == 4  # header + 3 rows

    def test_all_outputs_together(self, tmp_path) -> None:
        fasta = str(tmp_path / "out.fa")
        jsn = str(tmp_path / "out.json")
        tsv = str(tmp_path / "out.tsv")
        main([
            "--count", "5", "--length", "30",
            "--fasta", fasta, "--json", jsn, "--tsv", tsv,
            "--quiet", "--seed", "42",
        ])
        assert os.path.exists(fasta)
        assert os.path.exists(jsn)
        assert os.path.exists(tsv)

    def test_reproducible_with_seed(self, tmp_path) -> None:
        f1 = str(tmp_path / "a.fa")
        f2 = str(tmp_path / "b.fa")
        main(["--count", "5", "--fasta", f1, "--quiet", "--seed", "7"])
        main(["--count", "5", "--fasta", f2, "--quiet", "--seed", "7"])
        assert open(f1).read() == open(f2).read()

    def test_invalid_count_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--count", "0"])
        assert exc.value.code != 0

    def test_invalid_length_exits_nonzero(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--length", "0"])
        assert exc.value.code != 0

    def test_stdout_summary_printed(self, capsys) -> None:
        main(["--count", "2", "--length", "20", "--seed", "1"])
        captured = capsys.readouterr()
        assert "Generated 2 oligos" in captured.out

    def test_prefix_used_in_names(self, tmp_path) -> None:
        path = str(tmp_path / "out.json")
        main([
            "--count", "2", "--length", "10",
            "--prefix", "myoligo",
            "--json", path,
            "--quiet", "--seed", "1",
        ])
        data = json.loads(open(path).read())
        assert all(item["name"].startswith("myoligo") for item in data)


# ---------------------------------------------------------------------------
# read_json
# ---------------------------------------------------------------------------


class TestReadJson:
    def _write_oligo_analyses(self, tmp_path) -> tuple[list[OligoAnalysis], str]:
        analyses = [
            analyse_oligo(DNA("ACGT"), name="o1"),
            analyse_oligo(DNA("TTTTAAAA"), name="o2"),
        ]
        path = str(tmp_path / "out.json")
        write_json(analyses, path)
        return analyses, path

    def test_returns_list(self, tmp_path) -> None:
        _, path = self._write_oligo_analyses(tmp_path)
        result = read_json(path)
        assert isinstance(result, list)

    def test_correct_count(self, tmp_path) -> None:
        analyses, path = self._write_oligo_analyses(tmp_path)
        result = read_json(path)
        assert len(result) == len(analyses)

    def test_returns_oligo_analysis_instances(self, tmp_path) -> None:
        _, path = self._write_oligo_analyses(tmp_path)
        result = read_json(path)
        for obj in result:
            assert isinstance(obj, OligoAnalysis)

    def test_sequence_preserved(self, tmp_path) -> None:
        analyses, path = self._write_oligo_analyses(tmp_path)
        result = read_json(path)
        assert result[0].sequence == analyses[0].sequence
        assert result[1].sequence == analyses[1].sequence

    def test_name_preserved(self, tmp_path) -> None:
        analyses, path = self._write_oligo_analyses(tmp_path)
        result = read_json(path)
        assert result[0].name == "o1"
        assert result[1].name == "o2"

    def test_gc_content_preserved(self, tmp_path) -> None:
        analyses, path = self._write_oligo_analyses(tmp_path)
        result = read_json(path)
        assert result[0].gc_content == pytest.approx(analyses[0].gc_content)

    def test_empty_json_returns_empty_list(self, tmp_path) -> None:
        path = str(tmp_path / "empty.json")
        with open(path, "w") as fh:
            fh.write("[]\n")
        assert read_json(path) == []

    def test_reads_structured_oligos(self, tmp_path) -> None:
        from OligoDesign.structured import StructuredOligo, generate_palindromic_motif

        oligos = [generate_palindromic_motif(rng=random.Random(i)) for i in range(3)]
        for i, o in enumerate(oligos):
            o.name = f"p{i}"
        path = str(tmp_path / "structured.json")
        write_json(oligos, path)
        result = read_json(path)
        assert len(result) == 3
        for obj in result:
            assert isinstance(obj, StructuredOligo)

    def test_structured_sequence_preserved(self, tmp_path) -> None:
        from OligoDesign.structured import generate_palindromic_motif

        oligos = [generate_palindromic_motif(rng=random.Random(0))]
        oligos[0].name = "p0"
        path = str(tmp_path / "s.json")
        write_json(oligos, path)
        result = read_json(path)
        assert result[0].sequence == oligos[0].sequence

    def test_structured_oligo_type_preserved(self, tmp_path) -> None:
        from OligoDesign.structured import generate_inverted_repeat

        oligos = [generate_inverted_repeat(rng=random.Random(0))]
        path = str(tmp_path / "ir.json")
        write_json(oligos, path)
        result = read_json(path)
        assert result[0].oligo_type == "inverted_repeat"

    def test_invalid_json_raises(self, tmp_path) -> None:
        path = str(tmp_path / "bad.json")
        with open(path, "w") as fh:
            fh.write("not json at all")
        with pytest.raises(json.JSONDecodeError):
            read_json(path)

    def test_non_array_json_raises(self, tmp_path) -> None:
        path = str(tmp_path / "obj.json")
        with open(path, "w") as fh:
            json.dump({"key": "value"}, fh)
        with pytest.raises(ValueError, match="Expected a JSON array"):
            read_json(path)

    def test_mixed_collection_round_trip(self, tmp_path) -> None:
        """read_json must correctly handle a JSON file with both OligoAnalysis
        and StructuredOligo records — each record is typed independently."""
        from OligoDesign.structured import StructuredOligo, generate_palindromic_motif

        analysis = analyse_oligo(DNA("ACGTACGT"), name="random_oligo")
        structured = generate_palindromic_motif(half_length=6, rng=random.Random(0))
        structured.name = "pal_oligo"

        # write_json accepts any WritableOligo, so a mixed list is valid output
        path = str(tmp_path / "mixed.json")
        write_json([analysis, structured], path)

        result = read_json(path)
        assert len(result) == 2
        assert isinstance(result[0], OligoAnalysis)
        assert isinstance(result[1], StructuredOligo)
        assert result[0].sequence == analysis.sequence
        assert result[1].sequence == structured.sequence


# ---------------------------------------------------------------------------
# sequence_logo
# ---------------------------------------------------------------------------


class TestSequenceLogo:
    """Tests for sequence_logo.

    These tests are skipped automatically when the optional viz dependencies
    (logomaker / matplotlib) are not installed.
    """

    @pytest.fixture(autouse=True)
    def _require_viz(self):
        pytest.importorskip("logomaker")
        pytest.importorskip("matplotlib")

    def _oligos(self) -> list[OligoAnalysis]:
        return [
            analyse_oligo(DNA("ACGTACGT"), name=f"o{i}") for i in range(5)
        ]

    def test_creates_png_file(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo.png")
        sequence_logo(self._oligos(), path)
        assert os.path.exists(path)

    def test_png_file_non_empty(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo.png")
        sequence_logo(self._oligos(), path)
        assert os.path.getsize(path) > 0

    def test_reads_from_json_path(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        json_path = str(tmp_path / "oligos.json")
        png_path = str(tmp_path / "logo.png")
        write_json(self._oligos(), json_path)
        sequence_logo(json_path, png_path)
        assert os.path.exists(png_path)

    def test_probability_logo_type(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo_prob.png")
        sequence_logo(self._oligos(), path, logo_type="probability")
        assert os.path.exists(path)

    def test_information_logo_type(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo_info.png")
        sequence_logo(self._oligos(), path, logo_type="information")
        assert os.path.exists(path)

    def test_invalid_logo_type_raises(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo.png")
        with pytest.raises(ValueError, match="logo_type"):
            sequence_logo(self._oligos(), path, logo_type="invalid")

    def test_empty_list_raises(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo.png")
        with pytest.raises(ValueError):
            sequence_logo([], path)

    def test_object_without_sequence_raises(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo

        path = str(tmp_path / "logo.png")
        with pytest.raises(ValueError, match="sequence"):
            sequence_logo([object()], path)

    def test_structured_oligos_source(self, tmp_path) -> None:
        from OligoDesign.sequence_logo import sequence_logo
        from OligoDesign.structured import generate_palindromic_motif

        oligos = [
            generate_palindromic_motif(half_length=8, rng=random.Random(i))
            for i in range(5)
        ]
        path = str(tmp_path / "logo_structured.png")
        sequence_logo(oligos, path)
        assert os.path.exists(path)


# ---------------------------------------------------------------------------
# Edge-case validation for has_tandem_repeat (new)
# ---------------------------------------------------------------------------


class TestHasTandemRepeatValidation:
    def test_min_unit_zero_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_unit"):
            has_tandem_repeat(DNA("ATATAT"), min_unit=0)

    def test_min_unit_negative_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_unit"):
            has_tandem_repeat(DNA("ATATAT"), min_unit=-1)

    def test_max_unit_less_than_min_unit_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="max_unit"):
            has_tandem_repeat(DNA("ATATAT"), min_unit=4, max_unit=2)

    def test_min_count_one_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_count"):
            has_tandem_repeat(DNA("ATATAT"), min_count=1)

    def test_min_count_zero_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_count"):
            has_tandem_repeat(DNA("ATATAT"), min_count=0)

    def test_min_count_two_is_valid(self) -> None:
        # ACGTACGT = 2 copies of ACGT; min_count=2 should be valid and return True
        assert has_tandem_repeat(DNA("ACGTACGT"), min_unit=4, max_unit=4, min_count=2) is True


# ---------------------------------------------------------------------------
# Edge-case validation for find_complementary_pairs (new)
# ---------------------------------------------------------------------------


class TestFindComplementaryPairsValidation:
    def test_min_overlap_zero_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_overlap"):
            find_complementary_pairs([DNA("ACGT")], ["o1"], min_overlap=0)

    def test_min_overlap_negative_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="min_overlap"):
            find_complementary_pairs([DNA("ACGT")], ["o1"], min_overlap=-5)

    def test_names_length_mismatch_raises_value_error(self) -> None:
        oligos = [DNA("ACGT"), DNA("TTTT")]
        with pytest.raises(ValueError, match="same length"):
            find_complementary_pairs(oligos, ["o1"])

    def test_duplicate_names_raises_value_error(self) -> None:
        oligos = [DNA("ACGT"), DNA("TTTT")]
        with pytest.raises(ValueError, match="unique"):
            find_complementary_pairs(oligos, ["o1", "o1"])
