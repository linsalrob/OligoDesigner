"""Random oligonucleotide generation and sequence analysis."""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass, field
from typing import Optional, Protocol

from .dna import DNA


# ---------------------------------------------------------------------------
# Shared output protocol
# ---------------------------------------------------------------------------


class WritableOligo(Protocol):
    """Protocol for oligo objects that can be written to FASTA, JSON, and TSV.

    Any class that exposes ``name``, ``sequence``, ``to_dict()``,
    ``to_tsv_row()``, and the static ``tsv_headers()`` method satisfies this
    protocol and can be passed to :func:`write_fasta`, :func:`write_json`,
    and :func:`write_tsv`.
    """

    name: str
    sequence: str

    def to_dict(self) -> dict: ...

    def to_tsv_row(self) -> list[str]: ...

    @staticmethod
    def tsv_headers() -> list[str]: ...

# ---------------------------------------------------------------------------
# Random oligo generation
# ---------------------------------------------------------------------------

_BASES = "ACGT"


def random_oligo(length: int = 40, rng: Optional[random.Random] = None) -> DNA:
    """Generate a random DNA oligonucleotide of the given length.

    Parameters
    ----------
    length:
        Number of bases to generate.  Default is 40.
    rng:
        Optional :class:`random.Random` instance for reproducibility.
        If ``None`` a new instance is created (non-seeded).

    Returns
    -------
    DNA
        A new :class:`~OligoDesign.dna.DNA` object of the requested length.

    Raises
    ------
    ValueError
        If *length* is negative.

    Examples
    --------
    >>> import random
    >>> oligo = random_oligo(10, rng=random.Random(0))
    >>> len(oligo)
    10
    """
    if length < 0:
        raise ValueError(f"length must be >= 0, got {length}")
    if rng is None:
        rng = random.Random()
    return DNA("".join(rng.choice(_BASES) for _ in range(length)))


# ---------------------------------------------------------------------------
# Tandem-repeat detection
# ---------------------------------------------------------------------------


def has_tandem_repeat(
    dna: DNA,
    min_unit: int = 2,
    max_unit: int = 4,
    min_count: int = 3,
) -> bool:
    """Return ``True`` if the sequence contains a tandem repeat.

    Searches for consecutive occurrences of any motif whose length is
    between *min_unit* and *max_unit* (inclusive).

    Parameters
    ----------
    dna:
        The DNA sequence to test.
    min_unit:
        Minimum repeat-unit length in bases.  Default is 2.
    max_unit:
        Maximum repeat-unit length in bases.  Default is 4.
    min_count:
        Minimum number of consecutive copies required.  Default is 3.

    Returns
    -------
    bool

    Examples
    --------
    >>> has_tandem_repeat(DNA("ATATAT"))  # AT repeated 3 times
    True
    >>> has_tandem_repeat(DNA("ACGTACGT"))  # only 2 repeats of ACGT
    False
    """
    seq = str(dna)
    for unit_len in range(min_unit, max_unit + 1):
        repeat_len = unit_len * min_count
        for i in range(len(seq) - repeat_len + 1):
            unit = seq[i : i + unit_len]
            if seq[i : i + repeat_len] == unit * min_count:
                return True
    return False


# ---------------------------------------------------------------------------
# Cross-complementarity
# ---------------------------------------------------------------------------


def _complementarity_score(seq_a: str, seq_b: str, min_overlap: int) -> bool:
    """Return True if *seq_a* and *seq_b* share >= *min_overlap* complementary bases."""
    if len(seq_a) < min_overlap or len(seq_b) < min_overlap:
        return False
    # Build the reverse complement of seq_b.  Any substring of seq_a that
    # appears as a substring of rc_b represents a complementary region.
    rc_b = seq_b[::-1].translate(str.maketrans("ACGT", "TGCA"))
    # Pre-build a set of all min_overlap-length windows from rc_b so each
    # lookup in seq_a is O(k) rather than O(M).
    rc_b_windows = {
        rc_b[j : j + min_overlap]
        for j in range(len(rc_b) - min_overlap + 1)
    }
    for i in range(len(seq_a) - min_overlap + 1):
        if seq_a[i : i + min_overlap] in rc_b_windows:
            return True
    return False


def find_complementary_pairs(
    oligos: list[DNA],
    names: list[str],
    min_overlap: int = 10,
) -> dict[str, list[str]]:
    """Identify pairs of oligos that share complementary sequence.

    Parameters
    ----------
    oligos:
        List of :class:`~OligoDesign.dna.DNA` objects.
    names:
        Names for each oligo (same order as *oligos*).
    min_overlap:
        Minimum length of complementary overlap to flag.  Default is 10.

    Returns
    -------
    dict[str, list[str]]
        Mapping of oligo name to list of names of complementary partners.
    """
    result: dict[str, list[str]] = {name: [] for name in names}
    seqs = [str(o) for o in oligos]
    for i in range(len(oligos)):
        for j in range(i + 1, len(oligos)):
            if _complementarity_score(seqs[i], seqs[j], min_overlap):
                result[names[i]].append(names[j])
                result[names[j]].append(names[i])
    return result


# ---------------------------------------------------------------------------
# Per-oligo analysis
# ---------------------------------------------------------------------------


@dataclass
class OligoAnalysis:
    """Analysis results for a single oligonucleotide.

    Attributes
    ----------
    name:
        Sequence identifier.
    sequence:
        The DNA sequence string (uppercase).
    length:
        Sequence length in bases.
    gc_content:
        Fraction of G+C bases (0.0 – 1.0).
    entropy:
        Shannon entropy of the sequence in bits (0.0 – 2.0).
    tm:
        Estimated melting temperature in °C (nearest-neighbor model,
        50 mM Na⁺, 250 nM oligo concentration).  ``None`` if the
        sequence is too short for a meaningful calculation (< 2 bases).
    base_composition:
        Per-base counts ``{"A": n, "C": n, "G": n, "T": n}``.
    longest_homopolymer:
        Length of the longest single-base run.
    has_homopolymer:
        ``True`` if any homopolymer run meets the minimum length threshold.
    is_low_complexity:
        ``True`` if any sliding window is dominated by a single base.
    is_palindrome:
        ``True`` if the sequence equals its own reverse complement.
    has_hairpin:
        ``True`` if a potential hairpin structure was detected.
    has_tandem_repeat:
        ``True`` if a tandem repeat of a short motif was detected.
    complementary_to:
        Names of other oligos in the set that are cross-complementary.
    """

    name: str
    sequence: str
    length: int
    gc_content: float
    entropy: float
    tm: float | None
    base_composition: dict[str, int]
    longest_homopolymer: int
    has_homopolymer: bool
    is_low_complexity: bool
    is_palindrome: bool
    has_hairpin: bool
    has_tandem_repeat: bool
    complementary_to: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Return a plain-Python dictionary representation."""
        return asdict(self)

    def to_tsv_row(self) -> list[str]:
        """Return values as a list of strings, ready for TSV writing."""
        bc = self.base_composition
        return [
            self.name,
            self.sequence,
            str(self.length),
            f"{self.gc_content:.4f}",
            f"{self.entropy:.4f}",
            f"{self.tm:.2f}" if self.tm is not None else "",
            str(bc.get("A", 0)),
            str(bc.get("C", 0)),
            str(bc.get("G", 0)),
            str(bc.get("T", 0)),
            str(self.longest_homopolymer),
            str(self.has_homopolymer),
            str(self.is_low_complexity),
            str(self.is_palindrome),
            str(self.has_hairpin),
            str(self.has_tandem_repeat),
            ",".join(self.complementary_to),
        ]

    @staticmethod
    def tsv_headers() -> list[str]:
        """Return the TSV column headers."""
        return [
            "name",
            "sequence",
            "length",
            "gc_content",
            "entropy",
            "tm",
            "count_A",
            "count_C",
            "count_G",
            "count_T",
            "longest_homopolymer",
            "has_homopolymer",
            "is_low_complexity",
            "is_palindrome",
            "has_hairpin",
            "has_tandem_repeat",
            "complementary_to",
        ]


def analyse_oligo(
    dna: DNA,
    name: str = "",
    *,
    min_stem: int = 4,
    min_loop: int = 3,
    max_loop: int = 8,
    min_hp_run: int = 4,
    lc_window: int = 10,
    lc_threshold: float = 0.7,
    min_repeat_unit: int = 2,
    max_repeat_unit: int = 4,
    min_repeat_count: int = 3,
) -> OligoAnalysis:
    """Analyse a single oligo and return an :class:`OligoAnalysis`.

    Parameters
    ----------
    dna:
        The sequence to analyse.
    name:
        Identifier for the oligo.
    min_stem / min_loop / max_loop:
        Hairpin detection parameters (see :meth:`DNA.has_hairpin`).
    min_hp_run:
        Minimum homopolymer run length for :meth:`DNA.has_homopolymer`.
    lc_window / lc_threshold:
        Low-complexity detection parameters (see :meth:`DNA.is_low_complexity`).
    min_repeat_unit / max_repeat_unit / min_repeat_count:
        Tandem repeat detection parameters (see :func:`has_tandem_repeat`).
    """
    return OligoAnalysis(
        name=name,
        sequence=str(dna),
        length=len(dna),
        gc_content=dna.gc_content(),
        entropy=dna.entropy(),
        tm=dna.melting_temperature(),
        base_composition=dna.base_composition(),
        longest_homopolymer=dna.longest_homopolymer(),
        has_homopolymer=dna.has_homopolymer(min_length=min_hp_run),
        is_low_complexity=dna.is_low_complexity(
            window=lc_window, threshold=lc_threshold
        ),
        is_palindrome=dna.is_self_complementary(),
        has_hairpin=dna.has_hairpin(
            min_stem=min_stem, min_loop=min_loop, max_loop=max_loop
        ),
        has_tandem_repeat=has_tandem_repeat(
            dna,
            min_unit=min_repeat_unit,
            max_unit=max_repeat_unit,
            min_count=min_repeat_count,
        ),
    )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def write_fasta(analyses: list[WritableOligo], path: str) -> None:
    """Write oligo sequences to *path* in FASTA format.

    Accepts any list of objects that satisfy :class:`WritableOligo`
    (e.g. :class:`OligoAnalysis` or
    :class:`~OligoDesign.structured.StructuredOligo`).

    Parameters
    ----------
    analyses:
        List of oligo objects with ``name`` and ``sequence`` attributes.
    path:
        Output file path.
    """
    with open(path, "w") as fh:
        for a in analyses:
            fh.write(f">{a.name}\n{a.sequence}\n")


def write_json(analyses: list[WritableOligo], path: str) -> None:
    """Write the full analysis data structure to *path* as JSON.

    Accepts any list of objects that satisfy :class:`WritableOligo`.

    Parameters
    ----------
    analyses:
        List of oligo objects with a ``to_dict()`` method.
    path:
        Output file path.
    """
    data = [a.to_dict() for a in analyses]
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, allow_nan=False)
        fh.write("\n")


def write_tsv(analyses: list[WritableOligo], path: str) -> None:
    """Write per-oligo analysis to *path* as a tab-separated file.

    Accepts any list of objects that satisfy :class:`WritableOligo`.
    Column headers are retrieved from ``type(analyses[0]).tsv_headers()``,
    so the correct header row is used automatically for each oligo type.

    Parameters
    ----------
    analyses:
        List of oligo objects with ``to_tsv_row()`` and a class-level
        ``tsv_headers()`` method.  All objects must share the same TSV schema
        (i.e. return the same headers); mixing types with different schemas
        raises :exc:`TypeError`.
    path:
        Output file path.

    Raises
    ------
    TypeError
        If *analyses* contains objects with different TSV schemas (e.g. a
        mix of :class:`OligoAnalysis` and
        :class:`~OligoDesign.structured.StructuredOligo`).
    """
    if not analyses:
        open(path, "w").close()
        return

    headers = type(analyses[0]).tsv_headers()
    for i, item in enumerate(analyses[1:], start=1):
        item_headers = type(item).tsv_headers()
        if item_headers != headers:
            raise TypeError(
                f"Mixed TSV schemas detected: item 0 uses {type(analyses[0]).__name__!r} "
                f"schema but item {i} uses {type(item).__name__!r} schema. "
                "All items passed to write_tsv must share the same TSV schema."
            )

    with open(path, "w") as fh:
        fh.write("\t".join(headers) + "\n")
        for a in analyses:
            fh.write("\t".join(a.to_tsv_row()) + "\n")
