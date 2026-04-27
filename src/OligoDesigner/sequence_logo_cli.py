"""Command-line interface for generating sequence logo PNG images.

Usage
-----
::

    generate-sequence-logo INPUT OUTPUT [options]

Run ``generate-sequence-logo --help`` for the full option list.

The input can be a JSON file written by ``generate-oligos`` /
``generate-structured-oligos``, or a plain FASTA file.  The output is a
PNG image.

Examples
--------
::

    # From a JSON file produced by generate-oligos
    generate-oligos --count 50 --seed 1 --json oligos.json
    generate-sequence-logo oligos.json logo.png

    # With a custom logo type and title
    generate-sequence-logo oligos.json logo.png --logo-type information \\
        --title "Random oligos (n=50)"

    # From a FASTA file
    generate-sequence-logo sequences.fasta logo.png --format fasta

    # High-resolution probability logo
    generate-sequence-logo oligos.json logo.png \\
        --logo-type probability --dpi 300 --width 14 --height 4
"""

from __future__ import annotations

import argparse
import os
import sys
from types import SimpleNamespace


def _read_fasta(path: str) -> list[SimpleNamespace]:
    """Read sequences from a FASTA file.

    Returns a list of objects with a ``sequence`` attribute so they are
    compatible with :func:`~OligoDesigner.sequence_logo.sequence_logo`.

    Parameters
    ----------
    path:
        Path to the FASTA file.

    Returns
    -------
    list
        List of :class:`types.SimpleNamespace` objects, each with a
        ``sequence`` attribute containing the uppercase nucleotide string.

    Raises
    ------
    OSError
        If the file cannot be opened.
    """
    sequences: list[SimpleNamespace] = []
    current: list[str] = []

    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if current:
                    sequences.append(SimpleNamespace(sequence="".join(current)))
                    current = []
            elif line.strip():
                current.append(line.strip().upper())

    if current:
        sequences.append(SimpleNamespace(sequence="".join(current)))

    return sequences


def _detect_format(path: str) -> str:
    """Return ``'fasta'`` or ``'json'`` based on the file extension.

    Parameters
    ----------
    path:
        File path to inspect.

    Returns
    -------
    str
        ``'fasta'`` for ``.fa``, ``.fasta``, ``.fna``, or ``.fas``
        extensions; ``'json'`` otherwise.
    """
    _, ext = os.path.splitext(path.lower())
    if ext in {".fa", ".fasta", ".fna", ".fas"}:
        return "fasta"
    return "json"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="generate-sequence-logo",
        description=(
            "Generate a sequence logo PNG from a collection of DNA "
            "oligonucleotides.  Input can be a JSON file produced by "
            "generate-oligos / generate-structured-oligos, or a FASTA file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Input file: JSON (from generate-oligos) or FASTA format.",
    )
    parser.add_argument(
        "output",
        metavar="OUTPUT",
        help="Output PNG file path.",
    )

    # ---- Logo style options ------------------------------------------------
    logo = parser.add_argument_group("logo")
    logo.add_argument(
        "--stack-order",
        choices=["value", "alphabetical"],
        default="value",
        metavar="ORDER",
        dest="stack_order",
        help=(
            "Stacking order of bases at each position: "
            "'value' places the largest value at the bottom (default), "
            "'alphabetical' keeps bases in A–T order from the top."
        ),
    )
    logo.add_argument(
        "--logo-type",
        choices=["counts", "probability", "information"],
        default="counts",
        metavar="TYPE",
        help=(
            "Logo type: 'counts' (raw nucleotide counts), "
            "'probability' (fraction of each base), or "
            "'information' (information-content bits)."
        ),
    )
    logo.add_argument(
        "--title",
        default="",
        metavar="TEXT",
        help="Title to display above the logo.",
    )
    logo.add_argument(
        "--color-scheme",
        default="classic",
        metavar="SCHEME",
        help=(
            "logomaker color scheme. Common choices: 'classic', "
            "'base_pairing', 'NajafabadiEtAl2017'."
        ),
    )

    # ---- Figure size / resolution options ----------------------------------
    fig = parser.add_argument_group("figure")
    fig.add_argument(
        "--width",
        type=float,
        default=10.0,
        metavar="W",
        help="Figure width in inches.",
    )
    fig.add_argument(
        "--height",
        type=float,
        default=3.0,
        metavar="H",
        help="Figure height in inches.",
    )
    fig.add_argument(
        "--dpi",
        type=int,
        default=150,
        metavar="N",
        help="Output resolution in dots per inch.",
    )

    # ---- Input format ------------------------------------------------------
    inp = parser.add_argument_group("input")
    inp.add_argument(
        "--format",
        choices=["auto", "json", "fasta"],
        default="auto",
        metavar="FMT",
        dest="input_format",
        help=(
            "Input file format: 'auto' detects from the file extension, "
            "'json', or 'fasta'."
        ),
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for the ``generate-sequence-logo`` command.

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

    # ---- Validate input ----------------------------------------------------
    if not os.path.isfile(args.input):
        parser.error(f"Input file not found: {args.input!r}")

    # ---- Validate output directory exists ----------------------------------
    output_dir = os.path.dirname(os.path.abspath(args.output))
    if not os.path.isdir(output_dir):
        parser.error(f"Output directory does not exist: {output_dir!r}")

    # ---- Determine input format --------------------------------------------
    fmt = args.input_format
    if fmt == "auto":
        fmt = _detect_format(args.input)

    # ---- Load sequences ----------------------------------------------------
    if fmt == "fasta":
        fasta_sequences = _read_fasta(args.input)
        if not fasta_sequences:
            parser.error(f"No sequences found in FASTA file: {args.input!r}")
        logo_source: str | list = fasta_sequences
    else:
        # Pass file path directly; sequence_logo() handles JSON loading.
        logo_source = args.input

    # ---- Generate the logo -------------------------------------------------
    try:
        from .sequence_logo import sequence_logo
    except ImportError as exc:  # pragma: no cover
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    try:
        sequence_logo(
            logo_source,
            args.output,
            logo_type=args.logo_type,
            title=args.title,
            figsize=(args.width, args.height),
            color_scheme=args.color_scheme,
            dpi=args.dpi,
            stack_order=args.stack_order,
        )
    except (ImportError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Sequence logo saved to {args.output!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
