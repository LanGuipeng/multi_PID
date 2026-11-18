"""Complete geometric decomposition with certified atom outputs."""

import time
from fractions import Fraction as F
from mar_certified.channels import validate_channel
from mar_certified.intervals import clear_caches, certify_expression
from mar_certified.reference import full_lattice
from .benchmark import distribution
from .common import write_json, decimal_output


@decimal_output
def run(output_dir):
    trials = []
    for rep in range(3):
        clear_caches()
        d = distribution()
        start = time.perf_counter()
        full = full_lattice(d, maximum_sources=4, max_seconds=60)
        stage = time.perf_counter()
        widths = []
        for e in full["atoms"].values():
            (v, n) = certify_expression(e, F("1e-12"))
            widths.append(v.width)
        trials.append(
            {
                "total_seconds": time.perf_counter() - start,
                "geometry_seconds": full["timings"]["geometry"],
                "inversion_seconds": full["timings"]["inversion"],
                "certification_seconds": time.perf_counter() - stage,
                "max_width": str(max(widths)),
            }
        )
    es = [F(1, 8), F(1, 4), F(1, 3), F(2, 5)]
    for mask in range(1, 16):
        coords = d.coordinates(mask)
        marginal = d.marginal((0,) + coords)
        for state, p in marginal.items():
            expected = F(1, 2)
            for i, x in zip(coords, state[1:]):
                expected *= es[i - 1] if x != state[0] else 1 - es[i - 1]
            assert p == expected
        for t in (0, 1):
            validate_channel(d.channel(mask, t))
    out = {
        "trials": trials,
        "median_seconds": sorted((x["total_seconds"] for x in trials))[1],
        "nodes": 166,
        "required_source_subsets": 15,
        "tolerance": "1e-12",
        "normalization_and_all_marginals_verified": True,
    }
    write_json(output_dir, "full_baseline_results.json", out)
    return out
