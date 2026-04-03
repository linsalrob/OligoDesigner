"""Command-line interface for generating random oligonucleotides.

Usage
-----
::

    generate-oligos [--count N] [--length L] [--seed S]
                    [--fasta FILE] [--json FILE] [--tsv FILE]

Run ``generate-oligos --help`` for the full option list.
"""

from __future__ import annotations

import argparse
import random
import sys

from .dna import DNA
from .oligo import (
    OligoAnalysis,
    analyse_oligo,
    find_complementary_pairs,
    random_oligo,
    write_fasta,
    write_json,
    write_tsv,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="generate-oligos",
        description=(
            "Generate random oligonucleotide sequences and analyse them for "
            "hairpins, palindromes, tandem repeats, homopolymers, "
            "low-complexity regions, and cross-complementarity."
        ),
    )

    # Generation options
    gen = parser.add_argument_group("generation")
    gen.add_argument(
        "--count",
        "-n",
        type=int,
        default=10,
        metavar="N",
        help="Number of oligos to generate (default: 10).",
    )
    gen.add_argument(
        "--length",
        "-l",
        type=int,
        default=40,
        metavar="L",
        help="Length of each oligo in bases (default: 40).",
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
        default="oligo",
        metavar="PREFIX",
        help="Name prefix for generated oligos (default: 'oligo').",
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
        help="Write full analysis data to FILE in JSON format.",
    )
    out.add_argument(
        "--tsv",
        metavar="FILE",
        default=None,
        help="Write per-oligo analysis to FILE as a tab-separated table.",
    )
    out.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress the default summary printed to stdout.",
    )

    return parser


def _print_summary(analyses: list[OligoAnalysis]) -> None:
    """Print a human-readable summary to stdout."""
    flagged = [a for a in analyses if (
        a.has_homopolymer
        or a.is_low_complexity
        or a.is_palindrome
        or a.has_hairpin
        or a.has_tandem_repeat
        or a.complementary_to
    )]

    print(f"Generated {len(analyses)} oligos ({len(flagged)} flagged)\n")
    header = (
        f"{'Name':<20} {'Length':>6} {'GC%':>6} {'Entropy':>8}  "
        f"{'Homopol':>7} {'LowCpx':>6} {'Palind':>6} "
        f"{'Hairpin':>7} {'TandRep':>7} {'XCompl':>6}"
    )
    print(header)
    print("-" * len(header))
    for a in analyses:
        xcompl = "yes" if a.complementary_to else "no"
        print(
            f"{a.name:<20} {a.length:>6} {a.gc_content * 100:>5.1f}% {a.entropy:>8.4f}  "
            f"{'yes' if a.has_homopolymer else 'no':>7} "
            f"{'yes' if a.is_low_complexity else 'no':>6} "
            f"{'yes' if a.is_palindrome else 'no':>6} "
            f"{'yes' if a.has_hairpin else 'no':>7} "
            f"{'yes' if a.has_tandem_repeat else 'no':>7} "
            f"{xcompl:>6}"
        )


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``generate-oligos`` command.

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
    if args.length < 1:
        parser.error("--length must be >= 1")

    rng = random.Random(args.seed)

    # Generate oligos
    width = len(str(args.count))
    names: list[str] = []
    oligos: list[DNA] = []
    for i in range(1, args.count + 1):
        name = f"{args.prefix}{i:0{width}}"
        names.append(name)
        oligos.append(random_oligo(length=args.length, rng=rng))

    # Per-oligo analysis (without cross-complementarity)
    analyses: list[OligoAnalysis] = [
        analyse_oligo(
            oligos[i],
            name=names[i],
            min_stem=args.min_stem,
            min_loop=args.min_loop,
            max_loop=args.max_loop,
            min_hp_run=args.min_hp_run,
        )
        for i in range(args.count)
    ]

    # Cross-complementarity
    pairs = find_complementary_pairs(oligos, names, min_overlap=args.min_overlap)
    for a in analyses:
        a.complementary_to = pairs.get(a.name, [])

    # Output
    if args.fasta:
        write_fasta(analyses, args.fasta)
    if args.json:
        write_json(analyses, args.json)
    if args.tsv:
        write_tsv(analyses, args.tsv)

    if not args.quiet:
        _print_summary(analyses)

    return 0


if __name__ == "__main__":
    sys.exit(main())
