# Installation

## Prerequisites

- Python 3.10 or higher
- `pip` (included with Python 3.4+)

## Install from Source

Clone the repository and install it using `pip`:

```bash
git clone https://github.com/linsalrob/OligoDesign.git
cd OligoDesign
pip install .
```

This installs the `OligoDesign` package and registers two command-line tools:
- `generate-oligos`
- `generate-structured-oligos`

## Development Install

To install in editable mode (so that changes to the source are reflected immediately without reinstalling):

```bash
git clone https://github.com/linsalrob/OligoDesign.git
cd OligoDesign
pip install -e ".[dev]"
```

The `[dev]` extra installs [pytest](https://docs.pytest.org/) so you can run the test suite.

## Verify the Installation

After installation, confirm that the command-line tools are available:

```bash
generate-oligos --help
generate-structured-oligos --help
```

To run the test suite:

```bash
pytest
```

## Virtual Environment (recommended)

It is good practice to install into a virtual environment to avoid conflicts with other packages:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```
