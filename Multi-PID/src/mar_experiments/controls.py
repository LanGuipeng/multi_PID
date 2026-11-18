"""Analytical controls and complete-lattice reconstruction checks."""

import itertools
import time
from fractions import Fraction as F
from mar_certified.cases import analytic_cases
from mar_certified.distributions import Distribution
from mar_certified.exact import evaluate
from mar_certified.lattice import antichains
from mar_certified.verification import verify
from .common import write_json, decimal_output


@decimal_output
def run(output_dir):
    cases = analytic_cases()
    selected = {k: cases[k] for k in ["copy", "xor"]}
    selected["four_copy"] = {(t, t, t, t, t): F(1, 2) for t in (0, 1)}
    selected["four_unique"] = {
        (t, t, a, b, c): F(1, 16) for (t, a, b, c) in itertools.product((0, 1), repeat=4)
    }
    for eta in [F(0), F(1, 8)]:
        tab = {}
        for xs in itertools.product((0, 1), repeat=4):
            for e in (0, 1):
                p = F(1, 16) * (eta if e else 1 - eta)
                if p:
                    tab[(sum(xs) % 2 ^ e,) + xs] = p
        selected["four_parity_" + str(eta)] = tab
    selected["pair_redundancy"] = {
        (t, u, u ^ t, v, v ^ t): F(1, 8) for (t, u, v) in itertools.product((0, 1), repeat=3)
    }
    selected["symmetric_bsc"] = {
        (t,) + xs: F(1, 2)
        * F(1, 4) ** sum((x != t for x in xs))
        * F(3, 4) ** sum((x == t for x in xs))
        for t in (0, 1)
        for xs in itertools.product((0, 1), repeat=4)
    }
    out = {}
    for name, tab in selected.items():
        d = Distribution(tab)
        start = time.perf_counter()
        v = verify(d, "1e-18", 4)
        assert v["status"] == "success"
        rows = []
        for a in antichains(d.n_sources, 4):
            r = evaluate(d, a, tol="1e-18")
            assert r["exact_zero"] or r["interval"].lower > 0
            if not r["exact_zero"]:
                rows.append({"alpha": a, "value": r["approximate"]})
        out[name] = {
            "status": v["status"],
            "atoms": v["atoms"],
            "seconds": time.perf_counter() - start,
            "max_MI_residual": str(
                max((q["absolute_residual_upper"] for q in v["reconstruction"]))
            ),
            "max_inversion_residual": str(
                max((q["absolute_residual_upper"] for q in v["independent_geometry"]))
            ),
            "nonzero_atoms": rows,
        }
    write_json(output_dir, "controls_results.json", out)
    return out
