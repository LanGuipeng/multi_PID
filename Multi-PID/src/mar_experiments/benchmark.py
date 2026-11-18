"""Cold-cache query benchmarks for four asymmetric binary symmetric channels."""

import itertools
import platform
import sys
import time
from fractions import Fraction as F
from mar_certified.distributions import Distribution
from mar_certified.exact import evaluate
from mar_certified.intervals import clear_caches
from mar_certified.lattice import boundary_sets, antichains
from mar_certified.verification import verify
from mar_certified.channels import affine_line
from .common import write_json, decimal_output


def distribution():
    es = [F(1, 8), F(1, 4), F(1, 3), F(2, 5)]
    tab = {}
    for t in (0, 1):
        for xs in itertools.product((0, 1), repeat=4):
            p = F(1, 2)
            for x, e in zip(xs, es):
                p *= e if x != t else 1 - e
            tab[(t,) + xs] = p
    return Distribution(tab)


def intersections(d, alpha, kind):
    c = boundary_sets(alpha, d.n_sources) if kind == "atom" else ()
    number = 0
    for t, p in d.target_probabilities().items():
        ch = [d.channel(a, t) for a in alpha + c]
        knots = sorted({F(0), F(1)} | {a / b for cols in ch for (a, b) in cols if a and b})
        for j, left in enumerate(knots):
            right = knots[j + 1] if j + 1 < len(knots) else None
            lines = [
                affine_line(cols, (left + right) / 2 if right is not None else left + 1)
                for cols in ch
            ] + [(F(0), F(0))]
            pts = set()
            for (a, b), (cc, e) in itertools.combinations(lines, 2):
                if b != e:
                    v = (cc - a) / (b - e)
                    if v > left and (right is None or v < right):
                        pts.add(v)
            number += len(pts)
    return number


def compact(r):
    return {
        "value": r["approximate"],
        "width": str(r["interval"].width),
        "L": str(r["interval"].lower),
        "U": str(r["interval"].upper),
        "timings": r["timings"],
        "metrics": r["metrics"],
        "series_terms": r["series_terms"],
        "exact_zero": r["exact_zero"],
    }


@decimal_output
def run(output_dir):
    results = {
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "processor": platform.processor(),
        },
        "distribution": {
            "mechanism": "T fair; independent E_i; X_i=T xor E_i",
            "errors": ["1/8", "1/4", "1/3", "2/5"],
            "support": 32,
        },
        "queries": [],
    }
    qs = [(1, 2, 4, 8), (1,), (3, 12), (3, 5, 9, 6, 10, 12), (15,), (1, 6)]
    for a in qs:
        a = tuple(sorted(a))
        for kind in ["cumulative", "atom"]:
            trials = []
            for rep in range(3):
                clear_caches()
                d = distribution()
                r = evaluate(d, a, kind=kind, tol="1e-12", max_seconds=60)
                trials.append(r["timings"]["total"])
            rr = compact(r)
            rr.update(
                alpha=a,
                kind=kind,
                boundary=boundary_sets(a, 4),
                required_channels=len(a) + (len(boundary_sets(a, 4)) if kind == "atom" else 0),
                trials_seconds=trials,
                median_seconds=sorted(trials)[1],
                intersections=intersections(d, a, kind),
            )
            results["queries"].append(rr)
    clear_caches()
    d = distribution()
    start = time.perf_counter()
    v = verify(d, tol="1e-18", maximum_sources=4)
    assert v["status"] == "success"
    elapsed = time.perf_counter() - start
    results["full_verification"] = {
        "seconds": elapsed,
        "atoms": v["atoms"],
        "status": v["status"],
        "max_MI_residual_bound": str(
            max((z["absolute_residual_upper"] for z in v["reconstruction"]))
        ),
        "max_geometry_inversion_residual_bound": str(
            max((z["absolute_residual_upper"] for z in v["independent_geometry"]))
        ),
    }
    positives = []
    zeros = []
    for a in antichains(4, 4):
        r = evaluate(d, a, tol="1e-18", max_seconds=60)
        assert r["exact_zero"] or r["interval"].lower > 0
        (zeros if r["exact_zero"] else positives).append({"alpha": a, "value": r["approximate"]})
    results["atoms"] = positives + zeros
    results["positive_atoms"] = len(positives)
    results["structural_zero_atoms"] = len(zeros)
    write_json(output_dir, "benchmark_results.json", results)
    return results
