# Reproducing the results

The repository separates the recorded manuscript results in `results/reference/` from new runs in `results/local/`. Both use the same result filenames, allowing the plotting command to consume either directory.

## Installation

Use Python 3.10 or newer in an isolated environment, then install from the repository root:

```bash
python -m pip install -e ".[experiments,plots,dev]"
```

The core query engine uses the Python standard library. The `experiments` extra adds mpmath for independent numerical checks. The `plots` extra adds NumPy and Matplotlib. The `dev` extra supplies test and packaging tools.

For a saved run, retain the Python and dependency versions alongside its JSON records. Different versions, operating systems, and machine loads can change timing and the refinement level completed before a time limit.

## Experiment commands

```bash
mar-reproduce --suite quick --output results/local
mar-reproduce --suite full --output results/local
```

The quick suite evaluates the switching and ternary examples, analytical controls, four-source query benchmarks, and the complete geometric baseline. The full suite also runs the small-atom and endpoint precision experiments. A completed experiment suite can include individual budget-limited queries; their status and actual achieved width are part of the result.

| Record | Contents |
|---|---|
| `examples_results.json` | Switching distribution, dependent ternary target, independent-noise control, complete atom tables and target contributions |
| `small_atoms_results.json` | Vanishing-atom family, float errors, certified intervals, endpoint comparisons, and support cases |
| `controls_results.json` | Copy, XOR, pure four-source controls, pair redundancy, and identical BSCs |
| `benchmark_results.json` | Query structures, required channels, exact segments, intersections, stage timings, and full verification |
| `full_baseline_results.json` | Complete geometric decomposition and certification timings |
| `precision_results.json` | Exact and endpoint runs across tolerances, histories, allowances, and midpoint-error bounds |
| `summary.json` | New-run suite status, experiment durations, and budget-limited query locations |

The generated files retain the individual query budgets. Keep the complete JSON, rather than only its rounded `approximate` value, when comparing implementations.

## Recorded environment and budgets

The committed reference calculations were made on November 17, 2026 using Python 3.12.14 on macOS 26.5.2 ARM64. Independent checks used mpmath 1.2.1 with 80 decimal digits for the switching and ternary examples and 100 digits for the small-atom family.

Certification uses exact integer and fraction arithmetic. The reported logarithm order `N` counts series terms, rather than floating mantissa bits. The endpoint accumulation precision `B` counts dyadic rounding bits.

| Experiment | Per-query resource settings |
|---|---|
| Switching and ternary exact queries | 90 s |
| Four-source exact query and full-lattice baselines | 60 s |
| Four-source endpoint precision study | 20 s; 32,768 cells per completed level; 20 levels; 256 logarithm terms |
| Small-atom endpoint study | 12 s; 16,384 cells per completed level; 12 levels; 128 logarithm terms |

These limits describe the recorded runs. The actual completed levels are stored in the output. For example, the four-source endpoint run requesting width `1e-3` completed level 10, with 20,480 cells and width `0.0018719503911484`; the next level would exceed its cell budget. The small-atom endpoint run requesting `1e-8` completed 9,216 cells and returned width `0.0044725013`. Both enclosures are valid but fail the requested width.

Three cold-cache serial trials were used for the exact-query comparison and the complete-decomposition baseline. Query totals in the manuscript table are medians; their detailed stage timings are from the third trial. Precision and small-atom rows are single runs. Shared-workstation wall times are measurements of these runs, not portable performance guarantees.

The complete geometric baseline has a recorded median of `0.653725` s and includes all 166 cumulative values, Möbius inversion, and certification of all atoms. The separate `1.321601` s complete-validation measurement also computes direct atoms and checks consistency; it measures a different workload.

## Checks performed by the experiments

1. Exact normalization of every input distribution and consistency of every required marginal channel.
2. Nonnegative direct envelopes and structural certificates for zero atoms.
3. Agreement with geometric MAR cumulative values followed by complete Möbius inversion.
4. Reconstruction of mutual information for every nonempty source subset.
5. Inclusion of independent high-precision values and, for endpoint comparisons, narrower rational exact-integral certificates.
6. The actual rational width criterion for every result reported as attaining its tolerance.

The switching and ternary expressions are compared exactly in a prime-log basis. The four-source benchmark and analytical controls have comparison residuals that simplify to the zero rational-log expression. The small-atom family also has an analytic entropy identity, checked using rational certificates and 100-digit numerical evaluation.

## Figures

```bash
mar-plot --results results/reference --output docs/figures
mar-plot --results results/local --output results/local/figures
```

Plotting requires `examples_results.json` and `precision_results.json`; a full run supplies both. The command writes:

- `envelope_switching.pdf` and `.png`: source payoff curves and their taskwise lower envelopes.
- `ternary_atoms.pdf` and `.png`: all 18 atoms, split by prior-weighted target contribution.
- `precision_cost.pdf` and `.png`: endpoint convergence and measured time versus requested width.

PDF output is vector graphics; PNG output uses 600 dpi. Text uses Times New Roman when available, otherwise STIX. Mathematical labels use STIX. Fonts can change the figure's appearance without changing its data.

## Interpreting the scope

The largest supplied joint table has 32 positive rows, and complete decompositions cover at most four sources. Costs for larger inputs depend on marginal support, the uncovered boundary, and rational branch geometry. The endpoint implementation follows the manuscript's global schedule; an intersection-aware partition or verified compact-support shortcut would be a different computational experiment.

The interval guarantee is numerical and conditional on the supplied probability table. For empirical counts, it does not include sampling uncertainty about a population distribution.
