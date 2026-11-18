# Release validation

Version 0.1.0 was validated locally on November 18, 2026. These checks concern the packaged implementation; the manuscript's original measurements remain in `results/reference/`.

## Environment

| Component | Version |
|---|---|
| Operating system | macOS 26.5.2, ARM64 |
| Python | 3.12.14 |
| mpmath | 1.2.1 |
| NumPy | 2.5.5 |
| Matplotlib | 3.11.2 |
| pytest | 9.2.0 |
| Ruff | 0.17.5 |
| setuptools | 84.2.1 |
| build | 1.6.1 |

## Installation and checks

Both a source distribution and a wheel were built. The wheel was installed into an isolated environment, and the complete experiment command was run outside the source directory. The final test run used the installed wheel and passed all **17 tests**. Ruff lint and formatting checks passed for all 28 Python files.

The tests cover analytical copy and XOR atoms, rational switching geometry, ternary full Möbius inversion, mutual-information reconstruction, a positive atom below binary64 cancellation error, endpoint inclusion and budget status, degenerate support, input validation, JSON certificates, and CLI exit codes. Large rational endpoints round-trip without changing the interpreter's global integer-string limit; malformed noninteger or nonfinite rational fields are rejected.

## Full experiment run

```bash
mar-reproduce --suite full --output results/local
```

| Experiment | Wall time (s) |
|---|---:|
| Exact examples | 1.245 |
| Analytical controls | 2.725 |
| Four-source query benchmark | 2.311 |
| Complete-decomposition baseline | 1.633 |
| Small atoms and support cases | 15.985 |
| Precision study | 28.493 |
| Full suite | 52.395 |

All six experiment modules completed. Two endpoint queries reached their cell budgets before the requested tolerance: the small-atom request at `1e-8` and the four-source request at `1e-3`. Their valid, wider enclosures and failure status are retained in the results. Suite completion means that the experiments finished and their checks passed; it does not convert those queries into successful certificates at the requested width.

The switching and two ternary systems reproduce all 40 original atom expressions and certificates. The four-source benchmark reproduces all 166 atom expressions and the 12 queried certificates. Timing is measured anew and depends on the execution environment.

## Figures and repository checks

```bash
mar-plot --results results/reference --output docs/figures
```

All three figure pairs were regenerated using the installed plotting command and visually inspected. The outputs are vector PDFs and 600 dpi PNGs. The README uses repository-relative image paths, so the figures render on GitHub without external hosting.

The committed numerical records have SHA-256 checksums. Source, tests, examples, and documentation were checked for Chinese text and local account paths. Build output, interpreter caches, virtual environments, and local experimental runs are excluded from the repository archive.

The GitHub Actions workflow configures Python 3.10, 3.12, and 3.13 on Linux and Windows. Those remote jobs will run after upload; the completed local validation reported here used Python 3.12 on macOS.
