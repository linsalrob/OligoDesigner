"""Command-line interface for generating structured oligonucleotides.

Usage
-----
::

    generate-structured-oligos [--type {palindrome,inverted_repeat,at_rich,all}]
                               [--count N] [--half-length L]
                               [--outer-arm-length L] [--inner-half-length L]
                               [--spacer-length S]
                               [--seed S] [--prefix PREFIX]
                               [--fasta FILE] [--json FILE] [--tsv FILE]
                               [--quiet]

Run ``generate-structured-oligos --help`` for the full option list.
"""

from __future__ import annotations

import argparse
import random
import sys

from .oligo import write_fasta, write_json, write_tsv
from .structured import (
    StructuredOligo,
    generate_at_rich_palindrome,
    generate_inverted_repeat,
    generate_palindromic_motif,
)

# Ordered list of available structured-oligo types
_ALL_TYPES = ["palindrome", "inverted_repeat", "at_rich"]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="generate-structured-oligos",
        description=(
            "Generate structured oligonucleotides: palindromic motifs, "
            "inverted repeats, and AT-rich palindromes.  Optionally write "
            "results to FASTA, JSON, and/or TSV files."
        ),
    )

    # Generation options
    gen = parser.add_argument_group("generation")
    gen.add_argument(
        "--type",
        "-t",
        choices=_ALL_TYPES + ["all"],
        default="all",
        metavar="TYPE",
        help=(
            "Type of structured oligo to generate: "
            "'palindrome', 'inverted_repeat', 'at_rich', or 'all' "
            "(generates equal numbers of each type).  Default: all."
        ),
    )
    gen.add_argument(
        "--count",
        "-n",
        type=int,
        default=5,
        metavar="N",
        help=(
            "Number of oligos to generate per type.  "
            "With --type all, 3× this number are produced.  Default: 5."
        ),
    )
    gen.add_argument(
        "--half-length",
        type=int,
        default=6,
        metavar="L",
        help=(
            "Arm half-length in bases for palindromic motifs and AT-rich "
            "palindromes (total = 2 × L + spacer).  Default: 6."
        ),
    )
    gen.add_argument(
        "--outer-arm-length",
        type=int,
        default=8,
        metavar="L",
        help="Outer arm length for inverted-repeat oligos.  Default: 8.",
    )
    gen.add_argument(
        "--inner-half-length",
        type=int,
        default=6,
        metavar="L",
        help="Inner core half-length for inverted-repeat oligos.  Default: 6.",
    )
    gen.add_argument(
        "--spacer-length",
        type=int,
        default=2,
        metavar="S",
        help=(
            "Spacer length between arms (0 or 2–6 bp).  "
            "For AT-rich palindromes the spacer uses 'N' bases.  Default: 2."
        ),
    )
    gen.add_argument(
        "--seed",
        type=int,
        default=None,
        metavar="S",
        help="Random seed for reproducibility (default: unset).",
    )
    gen.add_argument(
        "--prefix",
        default="soligo",
        metavar="PREFIX",
        help="Name prefix for generated oligos (default: 'soligo').",
    )

    # Output options
    out = parser.add_argument_group("output")
    out.add_argument(
        "--fasta",
        metavar="FILE",
        default=None,
        help="Write oligo sequences to FILE in FASTA format.",
    )
    out.add_argument(
        "--json",
        metavar="FILE",
        default=None,
        help="Write full structured oligo data to FILE in JSON format.",
    )
    out.add_argument(
        "--tsv",
        metavar="FILE",
        default=None,
        help=(
            "Write per-oligo data to FILE as a tab-separated table "
            "(includes hairpin, palindrome, GC content, and more)."
        ),
    )
    out.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress the default summary printed to stdout.",
    )

    return parser


def _print_summary(oligos: list[StructuredOligo]) -> None:
    """Print a human-readable summary table to stdout."""
    header = (
        f"{'Name':<20} {'Length':>6} {'Type':<20} {'GC%':>6}  "
        f"{'Palind':>6} {'Hairpin':>7}"
    )
    print(f"Generated {len(oligos)} structured oligos\n")
    print(header)
    print("-" * len(header))
    for oligo in oligos:
        print(
            f"{oligo.name:<20} {oligo.length:>6} {oligo.oligo_type:<20} "
            f"{oligo.gc_content * 100:>5.1f}%  "
            f"{'yes' if oligo.is_palindrome else 'no':>6} "
            f"{'yes' if oligo.has_hairpin else 'no':>7}"
        )


def _generate_batch(
    oligo_type: str,
    count: int,
    args: argparse.Namespace,
    rng: random.Random,
    start_index: int,
    width: int,
) -> list[StructuredOligo]:
    """Generate *count* oligos of the given *oligo_type*."""
    oligos: list[StructuredOligo] = []
    for i in range(count):
        name = f"{args.prefix}{start_index + i:0{width}}"
        if oligo_type == "palindrome":
            oligo = generate_palindromic_motif(
                half_length=args.half_length,
                spacer_length=args.spacer_length,
                rng=rng,
            )
        elif oligo_type == "inverted_repeat":
            oligo = generate_inverted_repeat(
                inner_half_length=args.inner_half_length,
                outer_arm_length=args.outer_arm_length,
                inner_spacer_length=args.spacer_length,
                rng=rng,
            )
        else:  # at_rich
            oligo = generate_at_rich_palindrome(
                half_length=args.half_length,
                spacer_length=args.spacer_length,
                rng=rng,
            )
        oligo.name = name
        oligos.append(oligo)
    return oligos


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``generate-structured-oligos`` command.

    Parameters
    ----------
    argv:
        Argument list; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Exit code (0 = success, non-zero = error).
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Validate arguments
    if args.count < 1:
        parser.error("--count must be >= 1")
    if args.half_length < 1:
        parser.error("--half-length must be >= 1")
    if args.outer_arm_length < 1:
        parser.error("--outer-arm-length must be >= 1")
    if args.inner_half_length < 1:
        parser.error("--inner-half-length must be >= 1")
    if args.spacer_length < 0:
        parser.error("--spacer-length must be >= 0")

    rng = random.Random(args.seed)

    # Determine which types to generate
    types_to_generate = _ALL_TYPES if args.type == "all" else [args.type]

    total = args.count * len(types_to_generate)
    width = len(str(total))

    all_oligos: list[StructuredOligo] = []
    idx = 1
    for oligo_type in types_to_generate:
        batch = _generate_batch(oligo_type, args.count, args, rng, idx, width)
        all_oligos.extend(batch)
        idx += args.count

    # Write outputs using the shared write_* functions from oligo.py
    if args.fasta:
        write_fasta(all_oligos, args.fasta)
    if args.json:
        write_json(all_oligos, args.json)
    if args.tsv:
        write_tsv(all_oligos, args.tsv)

    if not args.quiet:
        _print_summary(all_oligos)

    return 0


if __name__ == "__main__":
    sys.exit(main())
