# Contributing to `ted-namer`

Thank you for your interest in contributing to **`ted-namer`**! We welcome bug reports, feature requests, documentation improvements, and pull requests from the structural biology and bioinformatics community.

---

## Reporting Bugs and Requesting Features

Please use the [GitHub Issues tracker](https://github.com/KarelBerka/TED-namer/issues) to report bugs or propose enhancements:
- **Bug Reports**: Include the `ted-namer` version (`0.1.0`), Python version, the exact TED domain ID(s) or CLI command used, and the expected vs. observed output.
- **Naming or CATH Classification Edge Cases**: If you encounter a TED domain where the assigned name, specificity tier, or CATH code could be improved, please include the TED ID (e.g., `AF-P04637-F1-model_v4_TED01`) so we can add it to our regression test suite.

---

## Development Setup

1. **Fork and clone the repository**:
   ```bash
   git clone https://github.com/KarelBerka/TED-namer.git
   cd TED-namer
   ```

2. **Install dependencies (including development/test extras)**:
   Using [`uv`](https://github.com/astral-sh/uv) (recommended):
   ```bash
   uv sync --extra dev
   ```
   Or using standard `pip` in a virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -e ".[dev]"
   ```

---

## Running the Test Suite & Accuracy Validation

Before submitting a Pull Request, ensure all unit and live integration tests pass:

```bash
uv run --extra dev pytest
```

If your changes modify the boundary overlap engine (`ted_namer/overlap.py`), the 4-tier specificity aggregator (`ted_namer/aggregator.py`), or the CATH classifier (`ted_namer/cath_classifier.py`), please also run the blind hold-out cross-validation suite on the 10,000-domain benchmark dataset:

```bash
uv run python scripts/evaluate_accuracy.py
```

Verify that:
- Non-`T` CATH precision remains $>90\%$ at every hierarchy depth ($C, C.A, C.A.T, C.A.T.H$), and
- Legacy PDB header artifacts and low-information strings remain suppressed.

---

## Submitting a Pull Request

1. Create a descriptive feature branch (`git checkout -b feat/your-feature` or `fix/issue-description`).
2. Add unit tests in `tests/` covering any new functionality or bug fixes.
3. Ensure all existing and new tests pass (`uv run --extra dev pytest`).
4. Open a Pull Request against `main` with a clear summary of your changes.

---

## Code of Conduct

Participation in this project is governed by our [Code of Conduct](CODE_OF_CONDUCT.md).
