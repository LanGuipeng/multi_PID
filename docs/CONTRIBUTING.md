# Contributing

Contributions should preserve the Shannon MAR functional and make their numerical guarantees explicit.

## Local checks

```bash
python -m pip install -e ".[experiments,plots,dev]"
python -m pytest tests
mar-reproduce --suite quick --output results/local
```

Run the full suite when a change affects endpoint refinement, logarithm certification, or small-atom arithmetic:

```bash
mar-reproduce --suite full --output results/local
```

## Changes to mathematical code

Use exact rational arithmetic for probability input, breakpoint ordering, intersections, and branch decisions. Add a focused regression case for a mathematical or interface change. Suitable checks include an analytic atom, geometric Möbius inversion, subset mutual-information reconstruction, or a verified interval inclusion.

Keep structural zero certificates distinct from intervals containing zero. A budget-limited query must retain its actual status and achieved width. High-precision floating evaluation is a comparison value, not a replacement for an interval proof.

Keep code, comments, docstrings, and documentation in English. Prefer short docstrings and comments that explain a non-obvious mathematical step.

## Results and figures

Write new experiments to `results/local/` or another untracked output directory. Preserve the committed `results/reference/` measurements unless the change intentionally updates the reference experiment. An update should include the environment, budgets, all result files, and regenerated figures. Distinguish timing changes from changes in mathematical output.

When reporting a bug, provide a minimal exact distribution, the antichain, query method, tolerance, resource limits, package versions, and the returned result. Avoid including local account paths or private data.

## Pull requests

Describe the mathematical or practical problem, the resulting behavior, and the checks performed. Keep unrelated refactoring separate from changes to the numerical algorithms.
