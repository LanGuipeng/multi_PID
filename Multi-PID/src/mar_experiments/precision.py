"""Precision-cost comparison with explicitly bounded endpoint resources."""

from fractions import Fraction as F
from mar_certified.exact import evaluate
from mar_certified.endpoint import evaluate as endpoint
from mar_certified.intervals import clear_caches
from mar_certified.lattice import boundary_sets
from .benchmark import distribution, compact
from .common import write_json, decimal_output


@decimal_output
def run(output_dir):
    a = (1, 6)
    out = {
        "alpha": a,
        "boundary": boundary_sets(a, 4),
        "distribution": "benchmark.distribution()",
        "endpoint_budget": {"seconds_per_query": 20, "max_cells": 32768, "max_level": 20},
        "runs": [],
    }
    truth = evaluate(distribution(), a, tol="1e-35")["interval"]
    for eps in ["1e-1", "3e-2", "1e-2", "3e-3", "1e-3", "1e-6", "1e-12"]:
        clear_caches()
        r = evaluate(distribution(), a, tol=eps, max_seconds=20)
        exact = compact(r)
        rr = {"epsilon": eps, "exact": exact}
        if F(eps) >= F(1, 1000):
            clear_caches()
            r = endpoint(distribution(), a, tol=eps, max_seconds=20, max_cells=32768, max_level=20)
            rr["endpoint"] = {k: v for (k, v) in r.items() if k not in ["expression", "history"]}
            rr["endpoint"]["interval"] = {
                "L": str(r["interval"].lower),
                "U": str(r["interval"].upper),
                "width": str(r["interval"].width),
            }
            rr["endpoint"]["contains_exact_certificate"] = (
                r["interval"].lower <= truth.lower and truth.upper <= r["interval"].upper
            )
            rr["endpoint"]["midpoint_error_upper"] = str(
                max(
                    abs((r["interval"].lower + r["interval"].upper) / 2 - truth.lower),
                    abs((r["interval"].lower + r["interval"].upper) / 2 - truth.upper),
                )
            )
            rr["history"] = [
                {
                    "level": h["level"],
                    "cells": h["grid_cells"],
                    "width": str(h["interval"].width),
                    "errors": h["error_budget"],
                }
                for h in r["history"]
            ]
            assert rr["endpoint"]["contains_exact_certificate"]
        out["runs"].append(rr)
        write_json(output_dir, "precision_results.json", out)
    certify_midpoints(out)
    write_json(output_dir, "precision_results.json", out)
    return out


def certify_midpoints(report):
    truth = evaluate(distribution(), tuple(report["alpha"]), tol="1e-45")["interval"]
    for row in report["runs"]:
        result = row["exact"]
        midpoint = (F(result["L"]) + F(result["U"])) / 2
        result["midpoint_error_upper"] = str(
            max(abs(midpoint - truth.lower), abs(midpoint - truth.upper))
        )
        result["contains_exact_certificate"] = (
            F(result["L"]) <= truth.lower <= truth.upper <= F(result["U"])
        )
        assert result["contains_exact_certificate"]
        assert F(result["width"]) <= F(row["epsilon"])
