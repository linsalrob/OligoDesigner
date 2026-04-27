"""Tests for OligoDesigner.sequence_logo and OligoDesigner.sequence_logo_cli."""

from __future__ import annotations

import os
import tempfile
from types import SimpleNamespace

import pytest

from OligoDesigner.sequence_logo import (
    _build_count_matrix,
    _BASES,
    _VALID_STACK_ORDERS,
    sequence_logo,
)
from OligoDesigner.sequence_logo_cli import _build_parser, main


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_oligos(*seqs: str) -> list[SimpleNamespace]:
    return [SimpleNamespace(sequence=s) for s in seqs]


# ---------------------------------------------------------------------------
# _build_count_matrix
# ---------------------------------------------------------------------------

def test_build_count_matrix_basic():
    df = _build_count_matrix(["ACGT"])
    assert list(df.columns) == ["A", "C", "G", "T"]
    assert list(df["A"]) == [1, 0, 0, 0]
    assert list(df["T"]) == [0, 0, 0, 1]


def test_build_count_matrix_empty_raises():
    with pytest.raises(ValueError, match="empty"):
        _build_count_matrix([])


def test_build_count_matrix_truncates():
    df = _build_count_matrix(["ACGT", "AC"])
    assert len(df) == 2  # truncated to shortest


# ---------------------------------------------------------------------------
# sequence_logo – invalid arguments
# ---------------------------------------------------------------------------

def test_invalid_logo_type():
    oligos = _make_oligos("ACGT", "ACGT")
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "logo.png")
        with pytest.raises(ValueError, match="logo_type"):
            sequence_logo(oligos, out, logo_type="bad")


def test_invalid_stack_order():
    oligos = _make_oligos("ACGT", "ACGT")
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "logo.png")
        with pytest.raises(ValueError, match="stack_order"):
            sequence_logo(oligos, out, stack_order="bad")


def test_empty_source_raises():
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "logo.png")
        with pytest.raises(ValueError, match="empty"):
            sequence_logo([], out)


# ---------------------------------------------------------------------------
# sequence_logo – stack_order="value" (default)
# ---------------------------------------------------------------------------

def test_stack_order_value_produces_png():
    oligos = _make_oligos("ACGT", "ACGT", "ACGT")
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "logo.png")
        sequence_logo(oligos, out, stack_order="value")
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0


# ---------------------------------------------------------------------------
# sequence_logo – stack_order="alphabetical"
# ---------------------------------------------------------------------------

def test_stack_order_alphabetical_produces_png():
    oligos = _make_oligos("ACGT", "ACGT", "ACGT")
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "logo.png")
        sequence_logo(oligos, out, stack_order="alphabetical")
        assert os.path.exists(out)
        assert os.path.getsize(out) > 0


def test_bases_alpha_order():
    # Confirms that the alphabetical column list is A, C, G, T — first column
    # drawn at the top by logomaker's "fixed" stack order.
    assert _BASES == ["A", "C", "G", "T"]


def test_valid_stack_orders_constant():
    assert _VALID_STACK_ORDERS == {"value", "alphabetical"}


# ---------------------------------------------------------------------------
# sequence_logo_cli – --stack-order argument
# ---------------------------------------------------------------------------

def test_cli_parser_stack_order_default():
    parser = _build_parser()
    args = parser.parse_args(["in.json", "out.png"])
    assert args.stack_order == "value"


def test_cli_parser_stack_order_alphabetical():
    parser = _build_parser()
    args = parser.parse_args(["in.json", "out.png", "--stack-order", "alphabetical"])
    assert args.stack_order == "alphabetical"


def test_cli_parser_stack_order_choices():
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["in.json", "out.png", "--stack-order", "invalid"])


def test_cli_stack_order_value(tmp_path):
    from OligoDesigner.oligo import analyse_oligo, write_json
    from OligoDesigner.dna import DNA

    seqs = ["ACGTACGT", "AACCGGTT", "TTGGCCAA", "GCTAGCTA", "TATACGCG"]
    oligos = [analyse_oligo(DNA(s), name=f"o{i}") for i, s in enumerate(seqs)]
    json_path = str(tmp_path / "oligos.json")
    out_path = str(tmp_path / "logo.png")
    write_json(oligos, json_path)

    rc = main([json_path, out_path, "--stack-order", "value"])
    assert rc == 0
    assert os.path.exists(out_path)


def test_cli_stack_order_alphabetical(tmp_path):
    from OligoDesigner.oligo import analyse_oligo, write_json
    from OligoDesigner.dna import DNA

    seqs = ["ACGTACGT", "AACCGGTT", "TTGGCCAA", "GCTAGCTA", "TATACGCG"]
    oligos = [analyse_oligo(DNA(s), name=f"o{i}") for i, s in enumerate(seqs)]
    json_path = str(tmp_path / "oligos.json")
    out_path = str(tmp_path / "logo.png")
    write_json(oligos, json_path)

    rc = main([json_path, out_path, "--stack-order", "alphabetical"])
    assert rc == 0
    assert os.path.exists(out_path)
