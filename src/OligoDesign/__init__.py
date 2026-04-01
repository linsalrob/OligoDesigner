"""OligoDesign – DNA oligonucleotide design library."""

from .dna import DNA
from .oligo import (
    OligoAnalysis,
    WritableOligo,
    analyse_oligo,
    find_complementary_pairs,
    has_tandem_repeat,
    random_oligo,
    write_fasta,
    write_json,
    write_tsv,
)
from .structured import (
    StructuredOligo,
    generate_at_rich_palindrome,
    generate_inverted_repeat,
    generate_palindromic_motif,
)

__all__ = [
    "DNA",
    # Random oligo generation and analysis
    "OligoAnalysis",
    "WritableOligo",
    "analyse_oligo",
    "find_complementary_pairs",
    "has_tandem_repeat",
    "random_oligo",
    # Structured oligo generation
    "StructuredOligo",
    "generate_palindromic_motif",
    "generate_inverted_repeat",
    "generate_at_rich_palindrome",
    # Shared output helpers (work with any WritableOligo)
    "write_fasta",
    "write_json",
    "write_tsv",
]
