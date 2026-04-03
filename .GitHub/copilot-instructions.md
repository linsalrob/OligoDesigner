# Copilot instructions for OligoDesigner

## Repository purpose
OligoDesigner is a Python repository for building a robust DNA oligonucleotide design library and, over time, an application for designing synthetic DNA oligos.

The current priority is to create a reliable core sequence API centred on a `DNA` object, along with analysis and validation methods needed for oligo design. Our initial biological target is the design and evaluation of 40 bp DNA oligonucleotides for synthesis.

This repository should be treated as a reusable scientific software package. Correctness, explicit assumptions, maintainability, and test coverage are more important than cleverness or premature optimisation.

## Main development goals
Focus development in this order:

1. Build a clean and dependable `DNA` object
2. Add common sequence manipulation methods
3. Add sequence quality and structure checks relevant to oligo design
4. Build oligo-specific candidate generation, filtering, and scoring
5. Keep the codebase suitable for later extension into CLI, browser, and WebAssembly/Emscripten use

## Core expectations
- Primary language: Python
- Test framework: `pytest`
- Package name: `OligoDesigner`
- Design style: library first, application second
- Keep biological logic separate from interface code
- Prefer explicit, deterministic behaviour
- Add tests for all public behaviour

## Repository conventions
Assume the repository will evolve into a normal Python package layout. Until a `pyproject.toml` exists, do not invent unnecessary packaging complexity, but write code that will fit naturally into a standard package structure.

Preferred layout:
- `OligoDesigner/src` for package code
- `tests/` for `pytest` tests
- `README.md` for overview and examples
- `.github/` for automation and Copilot instructions

If the actual repository structure changes, follow the established repository layout rather than creating a second structure in parallel.

## DNA object guidance
The `DNA` object is the foundation of this repository and should remain simple, predictable, and biologically sensible.

### Required behaviour
The `DNA` object should:
- store DNA sequences in uppercase
- validate input sequences
- raise clear errors for invalid characters unless ambiguity is intentionally supported
- support `len()`, string conversion, repr, equality, iteration, and slicing
- support reverse, complement, and reverse-complement operations
- support subsequence extraction
- support GC content and base composition
- support motif finding
- support checks for homopolymers and low-complexity regions
- support hairpin and self-complementarity analysis, initially using simple transparent heuristics

### Indexing rules
Indexing behaviour must be explicit and unambiguous.

- Python-native indexing must remain standard 0-based
- If 1-based indexing is supported, it must be provided through explicitly named methods such as:
  - `get_base_1()`
  - `slice_1()`
- Do not overload normal Python indexing with mixed or ambiguous semantics
- Do not hide coordinate transformations

## Oligo design guidance
The main applied goal is the design of synthesizable 40 bp DNA oligos.

Code for oligo design should support:
- generation of candidate oligos from longer templates
- configurable constraints
- transparent filtering of poor candidates
- simple, interpretable scoring systems

Important sequence properties to evaluate include:
- GC content extremes
- long homopolymers
- low-complexity sequence
- self-complementarity
- inverted repeats
- likely hairpins
- repetitive elements
- dimerisation risk

Forty base pairs is the current design target, but do not hard-code the entire repository so that only 40 bp oligos are ever possible.

## Hairpin and structure analysis
Structural analysis should begin with simple and explainable heuristics rather than overly complicated thermodynamic models.

Early implementations may include:
- identifying inverted repeats
- configurable minimum stem length
- configurable loop length
- a simple heuristic score for hairpin likelihood

Design these APIs so more advanced thermodynamic or folding models can be added later without breaking the higher-level interfaces.

## Coding standards
When generating or modifying code in this repository:

- use type hints for new public code
- write docstrings for public classes and functions
- prefer small focused methods over large monolithic classes
- raise informative exceptions
- avoid hidden global state
- minimise dependencies
- preserve deterministic behaviour
- keep core sequence logic independent from UI, web, or operating-system-specific code

Prefer readable code over compact clever code.

## Testing expectations
`pytest` is the required testing framework.

Every change to public behaviour should include tests. The `DNA` object in particular should have thorough test coverage.

Add tests for:
- valid and invalid sequence construction
- uppercase normalisation
- empty and short sequences where relevant
- indexing and slicing behaviour
- reverse complement correctness
- palindromic sequences
- GC content calculations
- motif search
- homopolymer detection
- low-complexity detection
- self-complementarity heuristics
- hairpin detection heuristics
- oligo candidate filtering and scoring

Where possible, use biologically meaningful test cases rather than only trivial examples.

## How to work in this repository
When making changes:
- inspect existing code and tests first
- preserve naming consistency and API clarity
- add tests in `tests/` alongside new functionality
- avoid broad refactors unless they clearly improve maintainability
- avoid introducing dependencies unless they are justified by a real need

If build or packaging infrastructure is missing, do not fabricate elaborate setup files unless explicitly asked. Keep the code ready for later packaging, but focus on library quality first.

## What to prioritise
Prioritise work in this order:
1. `DNA` class correctness
2. test coverage
3. core sequence analysis methods
4. oligo filtering and scoring
5. future extensibility

## What to avoid
- Do not introduce ambiguous indexing semantics
- Do not mix core sequence logic with interface code
- Do not hard-code assumptions that all future oligos must always be 40 bp
- Do not add thermodynamic complexity before the heuristic layer is working and tested
- Do not add unnecessary dependencies
- Do not write code without tests

## Copilot-specific guidance
When suggesting code for this repository:
- prefer complete, runnable, tested implementations
- include `pytest` tests with new behaviour
- preserve backwards compatibility where practical
- keep public APIs explicit and biologically interpretable
- suggest extension points when useful, but keep initial implementations simple
- assume this code may later be compiled or adapted for browser-compatible execution, so keep the core logic portable

The long-term goal is to build a trusted Python toolkit for DNA sequence manipulation and oligo design, starting from a strong `DNA` abstraction and growing into a broader oligo design platform.
