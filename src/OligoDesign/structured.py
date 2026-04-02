"""Structured oligonucleotide generation.

Provides functions to create three types of structured oligonucleotides:

- **Palindromic motifs**: two sequence halves that are reverse complements of
  each other, with an optional 2–6 bp spacer in between
  (e.g. ``ATGCGA TCGCAT`` or ``CGTACG CGTACG``).
- **Inverted repeats**: a palindromic inner core flanked on each side by an
  arm and its reverse complement
  (e.g. ``GCTAGTAC ATGCGA TT TCGCAT CGATCGTA``).
- **AT-rich inverted palindromes**: AT-only arms with an ``N`` spacer, common
  in phage and bacterial genomes
  (e.g. ``ATATTA NN TAATAT``).

All generator functions return a :class:`StructuredOligo` dataclass that
carries the full sequence, individual component parts, and pre-computed
structural analysis fields.

Use the shared output helpers :func:`~OligoDesign.oligo.write_fasta`,
:func:`~OligoDesign.oligo.write_json`, and
:func:`~OligoDesign.oligo.write_tsv` from :mod:`OligoDesign.oligo` to write
results to FASTA, JSON, or TSV files — the same functions used for random
oligos work here because :class:`StructuredOligo` satisfies the
:class:`~OligoDesign.oligo.WritableOligo` protocol.

Usage
-----
::

    import random
    from OligoDesign.structured import (
        generate_palindromic_motif,
        generate_inverted_repeat,
        generate_at_rich_palindrome,
    )
    from OligoDesign.oligo import write_fasta, write_json, write_tsv

    rng = random.Random(42)
    oligos = [
        generate_palindromic_motif(half_length=6, rng=rng),
        generate_inverted_repeat(rng=rng),
        generate_at_rich_palindrome(spacer_length=2, rng=rng),
    ]
    write_fasta(oligos, "structured.fa")
    write_json(oligos, "structured.json")
    write_tsv(oligos, "structured.tsv")
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Optional

from .dna import DNA, _COMPLEMENT_TABLE

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_ALL_BASES = "ACGT"
_AT_BASES = "AT"

# Recommended spacer lengths for structured oligos
SPACER_LENGTHS: tuple[int, ...] = (0, 2, 3, 4, 5, 6)

# Complement table used by the hairpin scanner (ACGT only; N has no mapping)
_HAIRPIN_COMPLEMENT_TABLE: dict[int, int] = str.maketrans("ACGT", "TGCA")
# Frozenset of unambiguous bases for stem validation
_ACGT_SET: frozenset[str] = frozenset("ACGT")


def _random_seq(length: int, bases: str, rng: random.Random) -> str:
    """Return a random sequence of *length* drawn from *bases*."""
    return "".join(rng.choice(bases) for _ in range(length))


def _rc(seq: str) -> str:
    """Return the reverse complement of an ACGT sequence string."""
    return seq[::-1].translate(_COMPLEMENT_TABLE)


# ---------------------------------------------------------------------------
# StructuredOligo dataclass
# ---------------------------------------------------------------------------


@dataclass
class StructuredOligo:
    """A structured oligonucleotide with defined internal symmetry.

    Satisfies the :class:`~OligoDesign.oligo.WritableOligo` protocol, so
    instances can be passed directly to :func:`~OligoDesign.oligo.write_fasta`,
    :func:`~OligoDesign.oligo.write_json`, and
    :func:`~OligoDesign.oligo.write_tsv`.

    Attributes
    ----------
    sequence:
        The full nucleotide sequence (uppercase; may contain ``N`` for
        AT-rich palindromes).
    oligo_type:
        One of ``'palindromic_motif'``, ``'inverted_repeat'``, or
        ``'at_rich_palindrome'``.
    left_arm:
        Left (5') arm sequence (ACGT only).
    right_arm:
        Right (3') arm sequence (ACGT only).  For correctly generated
        oligos this equals ``RC(left_arm)``.
    spacer:
        Spacer sequence between the two arms (may be ``'N' * n`` for
        AT-rich palindromes, or an ACGT sequence for other types).
    inner_left:
        Inner-left sequence for inverted repeats; empty string for other
        types.
    inner_right:
        Inner-right sequence for inverted repeats; empty string for other
        types.
    name:
        Optional identifier string.
    """

    sequence: str
    oligo_type: str
    left_arm: str
    right_arm: str
    spacer: str
    inner_left: str
    inner_right: str
    name: str = ""

    # ------------------------------------------------------------------
    # Computed properties
    # ------------------------------------------------------------------

    @property
    def length(self) -> int:
        """Total sequence length in bases."""
        return len(self.sequence)

    @property
    def is_palindrome(self) -> bool:
        """``True`` if the outer arms are reverse complements of each other.

        This is always ``True`` for correctly generated structured oligos.
        The spacer is not considered in this check.
        """
        return self.left_arm == _rc(self.right_arm)

    @property
    def inner_is_palindrome(self) -> bool:
        """``True`` if the inner arms (inverted-repeat type) are reverse
        complements of each other.  Always ``True`` for correctly generated
        inverted-repeat oligos; ``False`` for other types."""
        if not self.inner_left:
            return False
        return self.inner_left == _rc(self.inner_right)

    @property
    def gc_content(self) -> float:
        """GC fraction computed from ACGT bases only (``N`` bases ignored)."""
        acgt = [b for b in self.sequence if b in _ALL_BASES]
        if not acgt:
            return 0.0
        return sum(1 for b in acgt if b in "GC") / len(acgt)

    @property
    def entropy(self) -> float:
        """Shannon entropy in bits, computed from ACGT bases only (``N`` bases ignored)."""
        acgt_only = "".join(b for b in self.sequence if b in _ALL_BASES)
        if not acgt_only:
            return 0.0
        return DNA(acgt_only).entropy()

    @property
    def has_hairpin(self) -> bool:
        """``True`` if the sequence contains a potential hairpin structure.

        Performs a local hairpin scan that preserves the original sequence
        positions and treats ``N`` spacer bases as non-pairing loop bases.
        Only A/T/C/G complement pairing is permitted in stems; any ``N``
        base in a candidate stem position disqualifies that stem candidate
        so it is not reported as a hairpin.  ``N`` bases are freely allowed
        in loop positions.

        Default parameters: ``min_stem=4``, ``min_loop=3``, ``max_loop=8``.

        Notes on ambiguity handling
        ---------------------------
        - ``N`` bases may occupy any loop position without restriction.
        - ``N`` in a stem position prevents a stem base-pair from forming,
          so a candidate stem that contains any ``N`` is rejected.
        - This differs from the old approach (which stripped ``N`` before
          scanning), because here the original positional context is
          preserved: ``N`` spacer bases contribute to the loop length
          without collapsing the sequence geometry.
        """
        seq = self.sequence
        n = len(seq)
        min_stem, min_loop, max_loop = 4, 3, 8
        for stem_len in range(min_stem, n // 2 + 1):
            for loop_len in range(min_loop, max_loop + 1):
                required = 2 * stem_len + loop_len
                if required > n:
                    break
                for i in range(n - required + 1):
                    stem5 = seq[i : i + stem_len]
                    stem3 = seq[i + stem_len + loop_len : i + required]
                    # Reject if either stem contains non-ACGT bases (e.g. N)
                    if not frozenset(stem5).issubset(_ACGT_SET):
                        continue
                    if not frozenset(stem3).issubset(_ACGT_SET):
                        continue
                    rc_stem5 = stem5[::-1].translate(_HAIRPIN_COMPLEMENT_TABLE)
                    if stem3 == rc_stem5:
                        return True
        return False

    @property
    def has_tandem_repeat(self) -> bool:
        """``True`` if the ACGT-only portion contains a tandem repeat."""
        from .oligo import has_tandem_repeat as _htr

        acgt_only = "".join(b for b in self.sequence if b in _ALL_BASES)
        if not acgt_only:
            return False
        return _htr(DNA(acgt_only))

    @property
    def tm(self) -> float | None:
        """Estimated melting temperature in °C.

        Computed from the ACGT-only bases of the sequence (``N`` spacer
        bases are ignored), using the nearest-neighbor model at 50 mM
        Na⁺ and 250 nM oligo concentration.  Returns ``None`` if the
        ACGT-only portion is shorter than 2 bases.
        """
        acgt_only = "".join(b for b in self.sequence if b in _ALL_BASES)
        return DNA(acgt_only).melting_temperature() if len(acgt_only) >= 2 else None

    # ------------------------------------------------------------------
    # WritableOligo protocol implementation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Return a plain-Python dictionary representation."""
        return {
            "name": self.name,
            "sequence": self.sequence,
            "length": self.length,
            "oligo_type": self.oligo_type,
            "left_arm": self.left_arm,
            "right_arm": self.right_arm,
            "spacer": self.spacer,
            "inner_left": self.inner_left,
            "inner_right": self.inner_right,
            "is_palindrome": self.is_palindrome,
            "inner_is_palindrome": self.inner_is_palindrome,
            "gc_content": self.gc_content,
            "entropy": self.entropy,
            "tm": self.tm,
            "has_hairpin": self.has_hairpin,
            "has_tandem_repeat": self.has_tandem_repeat,
        }

    def to_tsv_row(self) -> list[str]:
        """Return values as a list of strings, ready for TSV writing."""
        return [
            self.name,
            self.sequence,
            str(self.length),
            self.oligo_type,
            self.left_arm,
            self.right_arm,
            self.spacer,
            self.inner_left,
            self.inner_right,
            str(self.is_palindrome),
            str(self.inner_is_palindrome),
            f"{self.gc_content:.4f}",
            f"{self.entropy:.4f}",
            f"{self.tm:.2f}" if self.tm is not None else "",
            str(self.has_hairpin),
            str(self.has_tandem_repeat),
        ]

    @staticmethod
    def tsv_headers() -> list[str]:
        """Return the TSV column headers."""
        return [
            "name",
            "sequence",
            "length",
            "oligo_type",
            "left_arm",
            "right_arm",
            "spacer",
            "inner_left",
            "inner_right",
            "is_palindrome",
            "inner_is_palindrome",
            "gc_content",
            "entropy",
            "tm",
            "has_hairpin",
            "has_tandem_repeat",
        ]


# ---------------------------------------------------------------------------
# Generator functions
# ---------------------------------------------------------------------------


def generate_palindromic_motif(
    half_length: int = 6,
    spacer_length: int = 0,
    *,
    rng: Optional[random.Random] = None,
) -> StructuredOligo:
    """Generate a palindromic motif oligonucleotide.

    Creates a sequence where the right half is the reverse complement of the
    left half, with an optional spacer in between.  The resulting oligo always
    has ``is_palindrome == True``.

    Parameters
    ----------
    half_length:
        Number of bases in each arm.  The total length without spacer is
        ``2 * half_length``.  Default is 6.
    spacer_length:
        Number of random ACGT bases between the two arms.  Should be one of
        0, 2, 3, 4, 5, or 6.  Default is 0.
    rng:
        Optional :class:`random.Random` instance for reproducibility.

    Returns
    -------
    StructuredOligo
        A new structured oligo of type ``'palindromic_motif'``.

    Raises
    ------
    ValueError
        If *half_length* < 1 or *spacer_length* is not in
        :data:`SPACER_LENGTHS` (``(0, 2, 3, 4, 5, 6)``).

    Examples
    --------
    >>> import random
    >>> oligo = generate_palindromic_motif(half_length=6, rng=random.Random(0))
    >>> oligo.is_palindrome
    True
    >>> len(oligo.sequence) == 12
    True
    """
    if half_length < 1:
        raise ValueError(f"half_length must be >= 1, got {half_length}")
    if spacer_length not in SPACER_LENGTHS:
        raise ValueError(
            f"spacer_length must be one of {SPACER_LENGTHS}, got {spacer_length}"
        )

    if rng is None:
        rng = random.Random()

    left_arm = _random_seq(half_length, _ALL_BASES, rng)
    right_arm = _rc(left_arm)
    spacer = _random_seq(spacer_length, _ALL_BASES, rng) if spacer_length > 0 else ""

    return StructuredOligo(
        sequence=left_arm + spacer + right_arm,
        oligo_type="palindromic_motif",
        left_arm=left_arm,
        right_arm=right_arm,
        spacer=spacer,
        inner_left="",
        inner_right="",
    )


def generate_inverted_repeat(
    inner_half_length: int = 6,
    outer_arm_length: int = 8,
    inner_spacer_length: int = 2,
    *,
    rng: Optional[random.Random] = None,
) -> StructuredOligo:
    """Generate an inverted repeat oligonucleotide.

    Produces a sequence with a palindromic inner core flanked by an outer arm
    and its reverse complement::

        [outer_arm] [inner_left] [inner_spacer] [inner_right] [RC(outer_arm)]

    Both the outer arms and the inner core satisfy ``is_palindrome`` /
    ``inner_is_palindrome``.

    Parameters
    ----------
    inner_half_length:
        Length of each half of the inner palindromic core.  Default is 6.
    outer_arm_length:
        Length of each outer flanking arm.  Default is 8.
    inner_spacer_length:
        Number of random ACGT bases in the spacer between the inner core
        halves.  Should be one of 0, 2, 3, 4, 5, or 6.  Default is 2.
    rng:
        Optional :class:`random.Random` instance for reproducibility.

    Returns
    -------
    StructuredOligo
        A new structured oligo of type ``'inverted_repeat'``.

    Raises
    ------
    ValueError
        If *inner_half_length* < 1, *outer_arm_length* < 1, or
        *inner_spacer_length* is not in :data:`SPACER_LENGTHS`
        (``(0, 2, 3, 4, 5, 6)``).

    Examples
    --------
    >>> import random
    >>> oligo = generate_inverted_repeat(rng=random.Random(0))
    >>> oligo.is_palindrome
    True
    >>> oligo.inner_is_palindrome
    True
    """
    if inner_half_length < 1:
        raise ValueError(f"inner_half_length must be >= 1, got {inner_half_length}")
    if outer_arm_length < 1:
        raise ValueError(f"outer_arm_length must be >= 1, got {outer_arm_length}")
    if inner_spacer_length not in SPACER_LENGTHS:
        raise ValueError(
            f"inner_spacer_length must be one of {SPACER_LENGTHS}, got {inner_spacer_length}"
        )

    if rng is None:
        rng = random.Random()

    outer_arm = _random_seq(outer_arm_length, _ALL_BASES, rng)
    rc_outer = _rc(outer_arm)

    inner_left = _random_seq(inner_half_length, _ALL_BASES, rng)
    inner_right = _rc(inner_left)
    inner_spacer = (
        _random_seq(inner_spacer_length, _ALL_BASES, rng)
        if inner_spacer_length > 0
        else ""
    )

    full_seq = outer_arm + inner_left + inner_spacer + inner_right + rc_outer

    return StructuredOligo(
        sequence=full_seq,
        oligo_type="inverted_repeat",
        left_arm=outer_arm,
        right_arm=rc_outer,
        spacer=inner_spacer,
        inner_left=inner_left,
        inner_right=inner_right,
    )


def generate_at_rich_palindrome(
    half_length: int = 6,
    spacer_length: int = 2,
    *,
    use_n_spacer: bool = True,
    rng: Optional[random.Random] = None,
) -> StructuredOligo:
    """Generate an AT-rich inverted palindrome.

    Creates an AT-only palindromic sequence with a spacer, common in phage
    and bacterial genomes::

        [AT_arm] [spacer] [RC(AT_arm)]

    When *use_n_spacer* is ``True`` (the default) the spacer consists of
    ``N`` ambiguity bases (e.g. ``ATATTA NN TAATAT``).  Set *use_n_spacer*
    to ``False`` to use random ACGT bases as the spacer instead.

    Parameters
    ----------
    half_length:
        Length of each AT-only arm.  Default is 6.
    spacer_length:
        Length of the spacer.  Must be one of :data:`SPACER_LENGTHS`
        (``(0, 2, 3, 4, 5, 6)``).  Default is 2.
    use_n_spacer:
        If ``True``, use ``'N'`` for all spacer positions.  Default is
        ``True``.
    rng:
        Optional :class:`random.Random` instance for reproducibility.

    Returns
    -------
    StructuredOligo
        A new structured oligo of type ``'at_rich_palindrome'``.

    Raises
    ------
    ValueError
        If *half_length* < 1 or *spacer_length* is not in
        :data:`SPACER_LENGTHS` (``(0, 2, 3, 4, 5, 6)``).

    Examples
    --------
    >>> import random
    >>> oligo = generate_at_rich_palindrome(half_length=6, spacer_length=2,
    ...                                    rng=random.Random(0))
    >>> oligo.is_palindrome
    True
    >>> set(oligo.left_arm).issubset({'A', 'T'})
    True
    """
    if half_length < 1:
        raise ValueError(f"half_length must be >= 1, got {half_length}")
    if spacer_length not in SPACER_LENGTHS:
        raise ValueError(
            f"spacer_length must be one of {SPACER_LENGTHS}, got {spacer_length}"
        )

    if rng is None:
        rng = random.Random()

    at_arm = _random_seq(half_length, _AT_BASES, rng)
    rc_arm = _rc(at_arm)
    spacer = (
        "N" * spacer_length
        if use_n_spacer
        else _random_seq(spacer_length, _ALL_BASES, rng)
    )

    return StructuredOligo(
        sequence=at_arm + spacer + rc_arm,
        oligo_type="at_rich_palindrome",
        left_arm=at_arm,
        right_arm=rc_arm,
        spacer=spacer,
        inner_left="",
        inner_right="",
    )
