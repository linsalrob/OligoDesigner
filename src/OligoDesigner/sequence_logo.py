"""Sequence logo generation from oligonucleotide collections.

Provides :func:`sequence_logo`, which reads a JSON file produced by
:func:`~OligoDesigner.oligo.write_json` (or accepts an already-loaded list
of oligo objects) and writes a PNG sequence logo to disk.

The logo is built using `logomaker <https://logomaker.readthedocs.io>`_ and
`matplotlib <https://matplotlib.org>`_.  Both packages must be installed::

    pip install logomaker matplotlib

Usage
-----
::

    from OligoDesigner.sequence_logo import sequence_logo

    # From a JSON file produced by write_json / generate-oligos
    sequence_logo("oligos.json", "logo.png")

    # From an already-loaded list of oligo objects
    from OligoDesigner.oligo import analyse_oligo, read_json
    from OligoDesigner.dna import DNA
    oligos = [analyse_oligo(DNA("ACGTACGT"), name="o1")]
    sequence_logo(oligos, "logo.png")
"""

from __future__ import annotations

import os
from typing import Union

_BASES = list("ACGT")


def _build_count_matrix(sequences: list[str]) -> "pandas.DataFrame":
    """Return a position × base count DataFrame.

    Sequences of different lengths are truncated to the length of the
    shortest sequence so that all positions are fully covered.

    Parameters
    ----------
    sequences:
        List of uppercase DNA sequence strings (ACGT only; ``N`` and other
        ambiguity bases are silently skipped at each position).

    Returns
    -------
    pandas.DataFrame
        Integer count matrix with columns ``['A', 'C', 'G', 'T']`` and one
        row per position.

    Raises
    ------
    ValueError
        If *sequences* is empty.
    """
    import pandas as pd

    if not sequences:
        raise ValueError("Cannot build a count matrix from an empty sequence list.")

    length = min(len(s) for s in sequences)
    counts: dict[str, list[int]] = {b: [0] * length for b in _BASES}

    for seq in sequences:
        for i, base in enumerate(seq[:length]):
            if base in counts:
                counts[base][i] += 1

    return pd.DataFrame(counts, index=range(length))


def sequence_logo(
    source: Union[str, list],
    output_path: str,
    *,
    logo_type: str = "counts",
    title: str = "",
    figsize: tuple[float, float] = (10.0, 3.0),
    color_scheme: str = "classic",
    dpi: int = 150,
) -> None:
    """Generate a PNG sequence logo from a collection of oligonucleotides.

    The input can be either the path to a JSON file previously written by
    :func:`~OligoDesigner.oligo.write_json` or an already-loaded list of oligo
    objects (anything with a ``sequence`` attribute, e.g.
    :class:`~OligoDesigner.oligo.OligoAnalysis` or
    :class:`~OligoDesigner.structured.StructuredOligo`).

    When sequences differ in length the logo is drawn over the region covered
    by *all* sequences (i.e. truncated to the shortest sequence length).

    Parameters
    ----------
    source:
        Either a file-system path (``str``) to a JSON file or a list of oligo
        objects that have a ``sequence`` attribute.
    output_path:
        Destination path for the PNG file.  Parent directories must already
        exist.
    logo_type:
        Type of logo to draw.  One of:

        * ``"counts"`` – raw nucleotide counts at each position (default).
        * ``"probability"`` – fraction of each base at each position.
        * ``"information"`` – information-content (bits) weighted heights.
    title:
        Optional title to display above the logo.
    figsize:
        ``(width, height)`` in inches.  Default is ``(10, 3)``.
    color_scheme:
        `logomaker` color scheme name.  Common choices: ``"classic"``
        (default), ``"base_pairing"``, ``"NajafabadiEtAl2017"``.
    dpi:
        Resolution of the output PNG in dots per inch.  Default is 150.

    Raises
    ------
    ImportError
        If ``logomaker`` or ``matplotlib`` is not installed.
    ValueError
        If *source* is an empty collection, or *logo_type* is not recognised.

    Examples
    --------
    >>> import os, tempfile, random
    >>> from OligoDesigner.oligo import analyse_oligo, write_json
    >>> from OligoDesigner.dna import DNA
    >>> from OligoDesigner.sequence_logo import sequence_logo
    >>> oligos = [analyse_oligo(DNA("ACGTACGT"), name=f"o{i}") for i in range(3)]
    >>> with tempfile.TemporaryDirectory() as d:
    ...     json_path = os.path.join(d, "oligos.json")
    ...     png_path  = os.path.join(d, "logo.png")
    ...     write_json(oligos, json_path)
    ...     sequence_logo(json_path, png_path)
    ...     os.path.exists(png_path)
    True
    """
    try:
        import logomaker
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "sequence_logo requires the 'logomaker' and 'matplotlib' packages. "
            "Install them with: pip install logomaker matplotlib"
        ) from exc

    _VALID_TYPES = {"counts", "probability", "information"}
    if logo_type not in _VALID_TYPES:
        raise ValueError(
            f"logo_type must be one of {sorted(_VALID_TYPES)!r}, got {logo_type!r}"
        )

    # ---- Load sequences ------------------------------------------------
    if isinstance(source, str):
        from .oligo import read_json

        objects = read_json(source)
    else:
        objects = list(source)

    if not objects:
        raise ValueError("Cannot draw a sequence logo from an empty collection.")

    sequences = []
    for i, obj in enumerate(objects):
        if not hasattr(obj, "sequence"):
            raise ValueError(
                f"Object at index {i} ({type(obj).__name__!r}) does not have a "
                "'sequence' attribute. Pass a list of OligoAnalysis or StructuredOligo "
                "objects, or a path to a JSON file written by write_json()."
            )
        sequences.append(obj.sequence.upper())

    # ---- Build count matrix --------------------------------------------
    count_df = _build_count_matrix(sequences)

    # ---- Convert to the requested logo type ----------------------------
    if logo_type == "counts":
        logo_df = count_df
    elif logo_type == "probability":
        row_sums = count_df.sum(axis=1)
        logo_df = count_df.div(row_sums, axis=0).fillna(0.0)
    else:  # information
        logo_df = logomaker.transform_matrix(
            count_df,
            from_type="counts",
            to_type="information",
        )

    # ---- Draw ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    logomaker.Logo(logo_df, ax=ax, color_scheme=color_scheme)

    ax.set_xlabel("Position")
    y_label = {"counts": "Count", "probability": "Probability", "information": "Bits"}
    ax.set_ylabel(y_label[logo_type])

    if title:
        ax.set_title(title)

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, format="png")
    plt.close(fig)
