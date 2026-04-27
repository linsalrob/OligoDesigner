"""Command-line interface for generating structured oligonucleotides.

Usage
-----
::

    generate-structured-oligos [--type {palindrome,inverted_repeat,at_rich,all}]
                               [--count N] [--half-length L]
                               [--outer-arm-length L] [--inner-half-length L]
                               [--spacer-length S]
                               [--five-prime-spacer SEQ] [--three-prime-spacer SEQ]
                               [--five-prime-random-length N] [--three-prime-random-length N]
                               [--seed S] [--prefix PREFIX]
                               [--min-stem N] [--min-loop N] [--max-loop N]
                               [--min-hp-run N] [--min-overlap N]
                               [--fasta FILE] [--json FILE] [--tsv FILE]
                               [--quiet]

Run ``generate-structured-oligos --help`` for the full option list.
"""

from __future__ import annotations

import argparse
import random
import sys

from .dna import DNA
from .oligo import find_complementary_pairs, remove_duplicate_sequences, write_fasta, write_json, write_tsv
from .structured import (
    SPACER_LENGTHS,
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

    # Flank / spacer options
    flank = parser.add_argument_group("flanks")
    flank.add_argument(
        "--five-prime-spacer",
        metavar="SEQ",
        default=None,
        help=(
            "ACGT sequence to prepend as a 5' flank to every oligo.  "
            "Mutually exclusive with --five-prime-random-length."
        ),
    )
    flank.add_argument(
        "--three-prime-spacer",
        metavar="SEQ",
        default=None,
        help=(
            "ACGT sequence to append as a 3' flank to every oligo.  "
            "Mutually exclusive with --three-prime-random-length."
        ),
    )
    flank.add_argument(
        "--five-prime-random-length",
        type=int,
        metavar="N",
        default=None,
        help=(
            "Generate a random ACGT sequence of length N as the 5' flank.  "
            "Mutually exclusive with --five-prime-spacer."
        ),
    )
    flank.add_argument(
        "--three-prime-random-length",
        type=int,
        metavar="N",
        default=None,
        help=(
            "Generate a random ACGT sequence of length N as the 3' flank.  "
            "Mutually exclusive with --three-prime-spacer."
        ),
    )

    # Analysis options
    ana = parser.add_argument_group("analysis")
    ana.add_argument(
        "--min-stem",
        type=int,
        default=4,
        metavar="N",
        help="Minimum stem length for hairpin detection (default: 4).",
    )
    ana.add_argument(
        "--min-loop",
        type=int,
        default=3,
        metavar="N",
        help="Minimum loop length for hairpin detection (default: 3).",
    )
    ana.add_argument(
        "--max-loop",
        type=int,
        default=8,
        metavar="N",
        help="Maximum loop length for hairpin detection (default: 8).",
    )
    ana.add_argument(
        "--min-hp-run",
        type=int,
        default=4,
        metavar="N",
        help="Minimum homopolymer run length to flag (default: 4).",
    )
    ana.add_argument(
        "--min-overlap",
        type=int,
        default=10,
        metavar="N",
        help="Minimum overlap for cross-complementarity detection (default: 10).",
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
    out.add_argument(
        "--deduplicate",
        action="store_true",
        help=(
            "Remove oligos with duplicate sequences before output, "
            "keeping only the first occurrence of each unique sequence."
        ),
    )

    return parser


def _resolve_flanks(args: argparse.Namespace, rng: random.Random) -> tuple[str, str]:
    """Return ``(five_prime_flank, three_prime_flank)`` strings from CLI args.

    Generates random ACGT sequences when the ``*_random_length`` arguments are
    provided; uses the literal user-supplied strings otherwise.  Returns empty
    strings when no flank option was given for that end.

    Parameters
    ----------
    args:
        Parsed argument namespace.
    rng:
        Random instance used when generating random flanks.

    Returns
    -------
    tuple[str, str]
        ``(five_prime_flank, three_prime_flank)`` – may be empty strings.
    """
    five_prime = ""
    three_prime = ""
    if args.five_prime_spacer is not None:
        five_prime = args.five_prime_spacer.upper()
    elif args.five_prime_random_length is not None:
        five_prime = "".join(
            rng.choice("ACGT") for _ in range(args.five_prime_random_length)
        )
    if args.three_prime_spacer is not None:
        three_prime = args.three_prime_spacer.upper()
    elif args.three_prime_random_length is not None:
        three_prime = "".join(
            rng.choice("ACGT") for _ in range(args.three_prime_random_length)
        )
    return five_prime, three_prime


def _print_summary(oligos: list[StructuredOligo]) -> None:
    """Print a human-readable summary table to stdout."""
    flagged = [o for o in oligos if (
        o.has_homopolymer
        or o.has_hairpin
        or o.has_tandem_repeat
        or o.complementary_to
    )]

    header = (
        f"{'Name':<20} {'Length':>6} {'Type':<20} {'GC%':>6} {'Entropy':>8}  "
        f"{'Palind':>6} {'Hairpin':>7} {'Homopol':>7} {'TandRep':>7} {'XCompl':>6}"
    )
    print(f"Generated {len(oligos)} structured oligos ({len(flagged)} flagged)\n")
    print(header)
    print("-" * len(header))
    for oligo in oligos:
        xcompl = "yes" if oligo.complementary_to else "no"
        print(
            f"{oligo.name:<20} {oligo.length:>6} {oligo.oligo_type:<20} "
            f"{oligo.gc_content * 100:>5.1f}% {oligo.entropy:>8.4f}  "
            f"{'yes' if oligo.is_palindrome else 'no':>6} "
            f"{'yes' if oligo.has_hairpin else 'no':>7} "
            f"{'yes' if oligo.has_homopolymer else 'no':>7} "
            f"{'yes' if oligo.has_tandem_repeat else 'no':>7} "
            f"{xcompl:>6}"
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
    if args.spacer_length not in SPACER_LENGTHS:
        parser.error(
            f"--spacer-length must be one of {list(SPACER_LENGTHS)}, "
            f"got {args.spacer_length}"
        )
    if args.min_stem < 1:
        parser.error("--min-stem must be >= 1")
    if args.min_loop < 1:
        parser.error("--min-loop must be >= 1")
    if args.max_loop < args.min_loop:
        parser.error("--max-loop must be >= --min-loop")
    if args.min_hp_run < 1:
        parser.error("--min-hp-run must be >= 1")
    if args.min_overlap < 1:
        parser.error("--min-overlap must be >= 1")

    # Flank validation
    if args.five_prime_spacer is not None and args.five_prime_random_length is not None:
        parser.error(
            "--five-prime-spacer and --five-prime-random-length are mutually exclusive"
        )
    if args.three_prime_spacer is not None and args.three_prime_random_length is not None:
        parser.error(
            "--three-prime-spacer and --three-prime-random-length are mutually exclusive"
        )
    if args.five_prime_spacer is not None and not set(
        args.five_prime_spacer.upper()
    ).issubset(set("ACGT")):
        parser.error("--five-prime-spacer must contain only A, C, G, T bases")
    if args.three_prime_spacer is not None and not set(
        args.three_prime_spacer.upper()
    ).issubset(set("ACGT")):
        parser.error("--three-prime-spacer must contain only A, C, G, T bases")
    if args.five_prime_random_length is not None and args.five_prime_random_length < 1:
        parser.error("--five-prime-random-length must be >= 1")
    if args.three_prime_random_length is not None and args.three_prime_random_length < 1:
        parser.error("--three-prime-random-length must be >= 1")

    rng = random.Random(args.seed)

    # Resolve flanks (before oligo generation to keep RNG state consistent)
    five_prime_flank, three_prime_flank = _resolve_flanks(args, rng)

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

    # Apply flanks if specified
    if five_prime_flank or three_prime_flank:
        for oligo in all_oligos:
            oligo.sequence = five_prime_flank + oligo.sequence + three_prime_flank

    # Remove duplicate sequences if requested
    if args.deduplicate:
        oligo_dnas = [DNA(o.sequence) for o in all_oligos]
        oligo_names = [o.name for o in all_oligos]
        _, unique_names, removed_names = remove_duplicate_sequences(oligo_dnas, oligo_names)
        unique_name_set = set(unique_names)
        all_oligos = [o for o in all_oligos if o.name in unique_name_set]
        if removed_names and not args.quiet:
            print(
                f"Removed {len(removed_names)} duplicate sequence(s): "
                + ", ".join(removed_names)
            )

    # Apply analysis parameters to each oligo
    for oligo in all_oligos:
        oligo.min_stem = args.min_stem
        oligo.min_loop = args.min_loop
        oligo.max_loop = args.max_loop
        oligo.min_hp_run = args.min_hp_run

    # Cross-complementarity (ACGT-only sequences used for matching)
    dna_seqs = [DNA("".join(b for b in o.sequence if b in "ACGT")) for o in all_oligos]
    names = [o.name for o in all_oligos]
    pairs = find_complementary_pairs(dna_seqs, names, min_overlap=args.min_overlap)
    for oligo in all_oligos:
        oligo.complementary_to = pairs.get(oligo.name, [])

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
