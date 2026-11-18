"""Command-line entry point for the paper's numerical experiments."""

import argparse
import importlib
import time
from pathlib import Path

from .common import write_json

QUICK = ("examples", "controls", "benchmark", "full_baseline")
SUITES = {
    "quick": QUICK,
    "full": QUICK + ("small_atoms", "precision"),
    "examples": ("examples",),
    "small-atoms": ("small_atoms",),
    "benchmark": ("benchmark", "full_baseline"),
    "precision": ("precision",),
    "controls": ("controls",),
}


def _budget_limited(data, path=""):
    rows = []
    if isinstance(data, dict):
        if data.get("status") == "timeout" and data.get("tolerance_met") is False:
            rows.append({"location": path, "reason": data.get("reason")})
        for key, value in data.items():
            rows.extend(_budget_limited(value, f"{path}/{key}"))
    elif isinstance(data, (list, tuple)):
        for index, value in enumerate(data):
            rows.extend(_budget_limited(value, f"{path}/{index}"))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mar-reproduce",
        description="Reproduce rational MAR experiments and write result JSON files.",
        epilog="Completed experiments retain budget-limited endpoint queries in their records.",
    )
    parser.add_argument("--suite", choices=tuple(SUITES), default="quick")
    parser.add_argument("--output", type=Path, default=Path("results/local"))
    arguments = parser.parse_args(argv)
    if not __debug__:
        parser.error("Run without -O: experiment validation uses assertions.")

    started = time.perf_counter()
    completed = []
    limited = []
    for name in SUITES[arguments.suite]:
        print(f"Running {name}...", flush=True)
        stage = time.perf_counter()
        try:
            module = importlib.import_module(f"mar_experiments.{name}")
            result = module.run(arguments.output)
        except ModuleNotFoundError as error:
            if error.name == "mpmath":
                parser.error('Install the experiment dependencies: pip install ".[experiments]"')
            raise
        completed.append({"experiment": name, "seconds": time.perf_counter() - stage})
        limited.extend(_budget_limited(result, name))
        print(f"Completed {name} in {completed[-1]['seconds']:.3f} s.", flush=True)
    write_json(
        arguments.output,
        "summary.json",
        {
            "suite": arguments.suite,
            "status": "completed",
            "experiments": completed,
            "budget_limited_queries": limited,
            "seconds": time.perf_counter() - started,
        },
    )
    print(f"Results written to {arguments.output}.", flush=True)
    if limited:
        print(
            f"{len(limited)} endpoint queries retained valid intervals without attaining their requested widths."
        )
