# OligoDesign

A Python library and command-line toolkit for designing DNA oligonucleotides for phage screens. OligoDesign generates random and structured oligos, analyses them for problematic secondary-structure features (hairpins, homopolymers, low-complexity regions, tandem repeats, and cross-complementarity), and writes results to FASTA, JSON, and TSV formats.

## Installation

See [INSTALLATION.md](INSTALLATION.md) for full instructions. In short:

```bash
git clone https://github.com/linsalrob/OligoDesign.git
cd OligoDesign
pip install .
```

Python 3.10 or higher is required.

## Command-line Tools

### `generate-oligos`

Generates and analyses random oligonucleotides.

```
generate-oligos [options]
```

#### Generation options

| Option | Default | Description |
|--------|---------|-------------|
| `--count N`, `-n N` | 10 | Number of oligos to generate |
| `--length L`, `-l L` | 40 | Length of each oligo in bases |
| `--seed S` | *(none)* | Integer seed for reproducible output |
| `--prefix PREFIX` | `oligo` | Name prefix for generated oligos |

#### Analysis options

| Option | Default | Description |
|--------|---------|-------------|
| `--min-stem N` | 4 | Minimum stem length for hairpin detection |
| `--min-loop N` | 3 | Minimum loop length for hairpin detection |
| `--max-loop N` | 8 | Maximum loop length for hairpin detection |
| `--min-hp-run N` | 4 | Minimum run length to call a homopolymer |
| `--min-overlap N` | 10 | Minimum overlap to call two oligos cross-complementary |

#### Output options

| Option | Description |
|--------|-------------|
| `--fasta FILE` | Write sequences to a FASTA file |
| `--json FILE` | Write full analysis results to a JSON file |
| `--tsv FILE` | Write analysis results to a tab-separated file |
| `--quiet`, `-q` | Suppress the summary table printed to stdout |

#### Examples

```bash
# Generate 20 oligos of length 35, reproducible, save as FASTA and TSV
generate-oligos --count 20 --length 35 --seed 42 --fasta output.fa --tsv output.tsv

# Stricter hairpin detection, save JSON only (no stdout output)
generate-oligos -n 5 -l 50 --min-stem 5 --max-loop 6 --quiet --json data.json
```

---

### `generate-structured-oligos`

Generates structured oligonucleotides with defined internal symmetry: palindromic motifs, inverted repeats, and AT-rich palindromes.

```
generate-structured-oligos [options]
```

#### Generation options

| Option | Default | Description |
|--------|---------|-------------|
| `--type TYPE`, `-t TYPE` | `all` | Type of oligo: `palindrome`, `inverted_repeat`, `at_rich`, or `all` |
| `--count N`, `-n N` | 5 | Number of oligos per type (with `--type all`, `3 × N` are produced) |
| `--half-length L` | 6 | Arm half-length (bp) for palindromes and AT-rich oligos |
| `--outer-arm-length L` | 8 | Outer arm length (bp) for inverted repeats |
| `--inner-half-length L` | 6 | Inner core half-length (bp) for inverted repeats |
| `--spacer-length S` | 2 | Spacer length: `0` or `2–6` bp |
| `--seed S` | *(none)* | Integer seed for reproducible output |
| `--prefix PREFIX` | `soligo` | Name prefix for generated oligos |

#### Output options

| Option | Description |
|--------|-------------|
| `--fasta FILE` | Write sequences to a FASTA file |
| `--json FILE` | Write full analysis results to a JSON file |
| `--tsv FILE` | Write analysis results to a tab-separated file |
| `--quiet`, `-q` | Suppress the summary table printed to stdout |

#### Oligo types

| Type | Structure | Always palindrome? |
|------|-----------|-------------------|
| `palindrome` | `[left_arm] [spacer] [RC(left_arm)]` | Yes |
| `inverted_repeat` | `[outer_arm] [inner_left] [spacer] [RC(inner_left)] [RC(outer_arm)]` | Yes |
| `at_rich` | AT-only arms with optional `N` spacer | Yes |

#### Examples

```bash
# One of each type, 3-bp spacer, save all formats
generate-structured-oligos --type all --count 3 --spacer-length 3 \
  --fasta output.fa --json output.json --tsv output.tsv

# Palindromic motifs only, no spacer, reproducible
generate-structured-oligos --type palindrome --count 10 --spacer-length 0 \
  --seed 123 --quiet --fasta palindromes.fa

# Inverted repeats with custom arm dimensions
generate-structured-oligos --type inverted_repeat -n 5 \
  --outer-arm-length 10 --inner-half-length 8 --spacer-length 2 \
  --prefix ir_ --fasta ir_oligos.fa
```

---

## Python API

OligoDesign is designed as a library first. The following public classes and functions are available after `pip install .`:

```python
from OligoDesign import (
    DNA,
    random_oligo,
    analyse_oligo,
    find_complementary_pairs,
    has_tandem_repeat,
    write_fasta, write_json, write_tsv,
    OligoAnalysis,
)
from OligoDesign.structured import (
    generate_palindromic_motif,
    generate_inverted_repeat,
    generate_at_rich_palindrome,
    StructuredOligo,
)
```

### `DNA` – core sequence object

```python
dna = DNA("ATCGATCGATCG")
dna.gc_content()          # float: 0.5
dna.entropy()             # float: Shannon entropy (bits)
dna.complement()          # DNA object
dna.reverse_complement()  # DNA object
dna.has_hairpin()         # bool (configurable stem/loop sizes)
dna.has_homopolymer()     # bool
dna.is_self_complementary()  # bool
dna.find_motif("ATCG")   # list[int] of 0-based positions
```

### Generating and analysing random oligos

```python
import random
rng = random.Random(42)

oligo = random_oligo(length=40, rng=rng)          # DNA object
analysis = analyse_oligo(oligo, name="oligo_001")  # OligoAnalysis dataclass

# Batch generation
oligos = [random_oligo(40, rng) for _ in range(100)]
names  = [f"oligo_{i:03d}" for i in range(100)]
analyses = [analyse_oligo(o, n) for o, n in zip(oligos, names)]

# Cross-complementarity check
pairs = find_complementary_pairs(oligos, names, min_overlap=10)

# Write output
write_fasta(analyses, "oligos.fa")
write_json(analyses, "oligos.json")
write_tsv(analyses, "oligos.tsv")
```

### Generating structured oligos

```python
rng = random.Random(0)

palindrome     = generate_palindromic_motif(half_length=6, spacer_length=2, rng=rng)
inverted       = generate_inverted_repeat(inner_half_length=6, outer_arm_length=8, rng=rng)
at_rich        = generate_at_rich_palindrome(half_length=6, spacer_length=2, rng=rng)

print(palindrome.sequence)       # e.g. "ATCGCG  CGCGAT"
print(palindrome.is_palindrome)  # True
```

## Output formats

| Format | Contents |
|--------|----------|
| FASTA  | Sequence name and nucleotide sequence |
| JSON   | Full per-oligo analysis (all fields from `OligoAnalysis` / `StructuredOligo`) |
| TSV    | Tab-separated table with one row per oligo; headers match JSON keys |

## License

MIT – see [LICENSE](LICENSE).
