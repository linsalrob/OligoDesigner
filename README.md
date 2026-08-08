# OligoDesigner

A Python library and command-line toolkit for designing DNA oligonucleotides for phage screens. OligoDesigner generates random and structured oligos, analyses them for problematic secondary-structure features (hairpins, homopolymers, low-complexity regions, tandem repeats, and cross-complementarity), and writes results to FASTA, JSON, and TSV formats.

## Browser application

The repository also includes a static, browser-only interface for all three tools. Open the [GitHub Pages site](https://linsalrob.github.io/OligoDesigner/) or serve the repository with any static file server and visit `/`. The application:

- generates random and structured oligos with the same options and analysis heuristics as the command-line tools;
- downloads selected FASTA, JSON, and TSV outputs locally;
- reads local FASTA or OligoDesigner JSON files and renders sequence logos as SVG; and
- has no server component, external runtime dependencies, analytics, or network requests. Input sequences never leave the browser.

To run the browser tests, use `npm run test:web` (Node.js 20 or newer; no install step is required).

## Installation

See [INSTALLATION.md](INSTALLATION.md) for full instructions. In short:

```bash
git clone https://github.com/linsalrob/OligoDesigner.git
cd OligoDesigner
pip install .
```

To also install the optional visualisation dependencies needed for sequence logo generation:

```bash
pip install ".[viz]"
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

#### Flank / spacer options

Spacers are added at the outer ends of every generated oligo.  For each end you can provide either an explicit sequence **or** a random-length option — not both.  If neither is given, no spacer is added.

| Option | Description |
|--------|-------------|
| `--five-prime-spacer SEQ` | Fixed ACGT sequence to prepend to every oligo (5' flank) |
| `--three-prime-spacer SEQ` | Fixed ACGT sequence to append to every oligo (3' flank) |
| `--five-prime-random-length N` | Generate a random ACGT sequence of length N as the 5' flank |
| `--three-prime-random-length N` | Generate a random ACGT sequence of length N as the 3' flank |

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

# Add fixed adapter sequences at both ends
generate-oligos --count 10 --length 40 \
  --five-prime-spacer AGATCGGAAGAGC \
  --three-prime-spacer CTCTTCCGATCT \
  --fasta flanked.fa

# Add random flanks of specified length (same length for each end but different sequences)
generate-oligos --count 10 --length 40 --seed 1 \
  --five-prime-random-length 8 --three-prime-random-length 8 \
  --fasta flanked.fa
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

#### Flank / spacer options

Spacers are added at the outer ends of every generated oligo.  For each end you can provide either an explicit sequence **or** a random-length option — not both.  If neither is given, no spacer is added.

| Option | Description |
|--------|-------------|
| `--five-prime-spacer SEQ` | Fixed ACGT sequence to prepend to every oligo (5' flank) |
| `--three-prime-spacer SEQ` | Fixed ACGT sequence to append to every oligo (3' flank) |
| `--five-prime-random-length N` | Generate a random ACGT sequence of length N as the 5' flank |
| `--three-prime-random-length N` | Generate a random ACGT sequence of length N as the 3' flank |

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

# Add fixed adapter sequences at both ends of every structured oligo
generate-structured-oligos --type palindrome --count 10 \
  --five-prime-spacer AGATCGGAAGAGC \
  --three-prime-spacer CTCTTCCGATCT \
  --fasta flanked.fa

# Add random flanks of specified length
generate-structured-oligos --type all --count 5 --seed 1 \
  --five-prime-random-length 8 --three-prime-random-length 8 \
  --fasta flanked.fa
```

---

### `generate-sequence-logo`

Generates a sequence logo PNG image from a collection of DNA oligonucleotides. Accepts a JSON file produced by `generate-oligos` or `generate-structured-oligos`, or a plain FASTA file.

```
generate-sequence-logo INPUT OUTPUT [options]
```

#### Positional arguments

| Argument | Description |
|----------|-------------|
| `INPUT` | Input file: JSON (from `generate-oligos` / `generate-structured-oligos`) or FASTA format |
| `OUTPUT` | Output PNG file path |

#### Logo options

| Option | Default | Description |
|--------|---------|-------------|
| `--logo-type TYPE` | `counts` | Logo style: `counts` (raw nucleotide counts), `probability` (fraction of each base), or `information` (information-content bits) |
| `--title TEXT` | *(none)* | Title to display above the logo |
| `--color-scheme SCHEME` | `classic` | [logomaker](https://logomaker.readthedocs.io) colour scheme (e.g. `classic`, `base_pairing`, `NajafabadiEtAl2017`) |

#### Figure options

| Option | Default | Description |
|--------|---------|-------------|
| `--width W` | `10.0` | Figure width in inches |
| `--height H` | `3.0` | Figure height in inches |
| `--dpi N` | `150` | Output resolution in dots per inch |

#### Input options

| Option | Default | Description |
|--------|---------|-------------|
| `--format FMT` | `auto` | Input format: `auto` (detected from file extension), `json`, or `fasta` |

#### Examples

```bash
# Generate oligos, save to JSON, then create a sequence logo
generate-oligos --count 50 --length 40 --seed 1 --json oligos.json
generate-sequence-logo oligos.json logo.png

# Information-content logo with a custom title
generate-sequence-logo oligos.json logo.png --logo-type information \
  --title "Random oligos (n=50)"

# From a FASTA file, high-resolution probability logo
generate-sequence-logo sequences.fasta logo.png --format fasta \
  --logo-type probability --dpi 300 --width 14 --height 4

# From structured oligos JSON
generate-structured-oligos --type palindrome --count 20 --seed 42 --json palindromes.json
generate-sequence-logo palindromes.json palindromes_logo.png --logo-type information
```

> **Note:** The `viz` extras must be installed: `pip install ".[viz]"`

---

## Python API

OligoDesigner is designed as a library first. The following public classes and functions are available after `pip install .`:

```python
from OligoDesigner import (
    DNA,
    random_oligo,
    analyse_oligo,
    find_complementary_pairs,
    has_tandem_repeat,
    write_fasta, write_json, write_tsv,
    read_json,
    OligoAnalysis,
    sequence_logo,
)
from OligoDesigner.structured import (
    generate_palindromic_motif,
    generate_inverted_repeat,
    generate_at_rich_palindrome,
    StructuredOligo,
)
```

### `DNA` – core sequence object

```python
dna = DNA("ATCGATCGATCG")
dna.gc_content()             # float: fraction of G+C bases (0.0–1.0)
dna.entropy()                # float: Shannon entropy in bits (0.0–2.0)
dna.melting_temperature()    # float | None: nearest-neighbor Tm in °C
dna.base_composition()       # dict[str, int]: {"A": n, "C": n, "G": n, "T": n}
dna.longest_homopolymer()    # int: length of longest single-base run
dna.complement()             # DNA object
dna.reverse_complement()     # DNA object
dna.has_hairpin()            # bool (configurable stem/loop sizes)
dna.has_homopolymer()        # bool
dna.is_low_complexity()      # bool (configurable window and threshold)
dna.is_self_complementary()  # bool
dna.find_motif("ATCG")      # list[int] of 0-based positions
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

# Read a previously written JSON file back into OligoAnalysis objects
analyses = read_json("oligos.json")
```

#### `OligoAnalysis` fields

| Field | Type | Description |
|-------|------|-------------|
| `name` | `str` | Sequence identifier |
| `sequence` | `str` | DNA sequence (uppercase) |
| `length` | `int` | Sequence length in bases |
| `gc_content` | `float` | Fraction of G+C bases (0.0–1.0) |
| `entropy` | `float` | Shannon entropy in bits (0.0–2.0) |
| `tm` | `float \| None` | Nearest-neighbor melting temperature in °C |
| `base_composition` | `dict` | Per-base counts `{"A": n, "C": n, "G": n, "T": n}` |
| `longest_homopolymer` | `int` | Length of the longest single-base run |
| `has_homopolymer` | `bool` | `True` if a homopolymer run meets the minimum length threshold |
| `is_low_complexity` | `bool` | `True` if any sliding window is dominated by a single base |
| `is_palindrome` | `bool` | `True` if the sequence equals its own reverse complement |
| `has_hairpin` | `bool` | `True` if a potential hairpin structure was detected |
| `has_tandem_repeat` | `bool` | `True` if a tandem repeat of a short motif was detected |
| `complementary_to` | `list[str]` | Names of cross-complementary oligos in the same set |

### Generating structured oligos

```python
rng = random.Random(0)

palindrome     = generate_palindromic_motif(half_length=6, spacer_length=2, rng=rng)
inverted       = generate_inverted_repeat(inner_half_length=6, outer_arm_length=8, rng=rng)
at_rich        = generate_at_rich_palindrome(half_length=6, spacer_length=2, rng=rng)

print(palindrome.sequence)       # e.g. "ATCGCG  CGCGAT"
print(palindrome.is_palindrome)  # True
```

## Sequence logo images

Sequence logos visualise the nucleotide frequency (or information content) at each position across a collection of oligos.  The `sequence_logo` function reads a JSON file produced by `write_json` (or accepts an already-loaded list of oligo objects) and saves a PNG image.

### Prerequisites

Install the visualisation extras if you haven't already:

```bash
pip install ".[viz]"
# or individually:
pip install logomaker matplotlib
```

### Python API

```python
from OligoDesigner import sequence_logo

# From a JSON file written by generate-oligos or write_json
sequence_logo("oligos.json", "logo.png")

# From an in-memory list of OligoAnalysis or StructuredOligo objects
sequence_logo(analyses, "logo.png")

# Customise the logo
sequence_logo(
    "oligos.json",
    "logo.png",
    logo_type="information",   # "counts" | "probability" | "information"
    title="My oligo set",
    figsize=(12, 4),           # width × height in inches
    color_scheme="classic",    # logomaker colour scheme
    dpi=200,
)
```

#### `sequence_logo` parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `source` | *(required)* | Path to a JSON file **or** a list of `OligoAnalysis` / `StructuredOligo` objects |
| `output_path` | *(required)* | Destination path for the PNG file |
| `logo_type` | `"counts"` | Logo style: `"counts"`, `"probability"`, or `"information"` (bits) |
| `title` | `""` | Optional title displayed above the logo |
| `figsize` | `(10.0, 3.0)` | Figure size in inches `(width, height)` |
| `color_scheme` | `"classic"` | [logomaker](https://logomaker.readthedocs.io) colour scheme name |
| `dpi` | `150` | Output resolution in dots per inch |

#### Logo types

| `logo_type` | Y-axis | Description |
|-------------|--------|-------------|
| `"counts"` | Count | Raw nucleotide counts at each position |
| `"probability"` | Probability | Fraction of each base at each position |
| `"information"` | Bits | Information-content weighted heights (max 2 bits per position) |

#### Example workflow

```bash
# Step 1 – generate oligos and save to JSON
generate-oligos --count 50 --length 40 --seed 1 --json oligos.json

# Step 2 – create sequence logo from the JSON file (Python)
python - <<'EOF'
from OligoDesigner import sequence_logo
sequence_logo("oligos.json", "logo.png", logo_type="information", title="Random oligos")
EOF
```

---



## Output formats

| Format | Contents |
|--------|----------|
| FASTA  | Sequence name and nucleotide sequence |
| JSON   | Full per-oligo analysis (all fields from `OligoAnalysis` / `StructuredOligo`) |
| TSV    | Tab-separated table with one row per oligo; headers match JSON keys |

## License

MIT – see [LICENSE](LICENSE).
