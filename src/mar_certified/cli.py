"""Command-line queries and complete-lattice verification."""

import argparse
import json
from pathlib import Path

from .api import query, serialize
from .distributions import Distribution


def _state_value(value):
    return tuple(_state_value(item) for item in value) if isinstance(value, list) else value


def load_input(path, cache=True):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    table = {}
    for row in data["records"]:
        state = tuple(_state_value(value) for value in row["state"])
        if state in table:
            raise ValueError("Duplicate states in input file")
        table[state] = row.get("mass", row.get("probability", row.get("count")))
    return Distribution(
        table, kind=data["kind"], float_policy=data.get("float_policy"), cache=cache
    )


def _nonnegative_float(value):
    result = float(value)
    if not result >= 0 or result == float("inf"):
        raise argparse.ArgumentTypeError("Must be a finite nonnegative number")
    return result


def _nonnegative_integer(value):
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError("Must be a nonnegative integer")
    return result


def _positive_integer(value):
    result = _nonnegative_integer(value)
    if not result:
        raise argparse.ArgumentTypeError("Must be a positive integer")
    return result


def parser():
    command = argparse.ArgumentParser(prog="mar-certified")
    sub = command.add_subparsers(dest="command", required=True)
    request = sub.add_parser("query", help="Evaluate a selected MAR cumulative value or atom")
    request.add_argument("input", help="JSON probability or count table")
    request.add_argument("--alpha", required=True, help="Source sets, e.g. '1|2,3'")
    request.add_argument("--kind", choices=["atom", "cumulative"], default="atom")
    request.add_argument(
        "--method", choices=["exact", "endpoint", "reference", "float"], default="exact"
    )
    request.add_argument("--tol", default="1e-8", help="Positive rational interval-width tolerance")
    request.add_argument("--algorithm", choices=["enumerate", "optimized"], default="enumerate")
    request.add_argument("--max-seconds", type=_nonnegative_float, default=20)
    request.add_argument("--max-cells", type=_nonnegative_integer, default=4096)
    request.add_argument("--max-level", type=_nonnegative_integer, default=20)
    request.add_argument("--max-terms", type=_positive_integer, default=256)
    request.add_argument("--maximum-sources", type=int, choices=range(1, 5), default=3)
    request.add_argument("--output", help="Write the JSON result to this path")
    verification = sub.add_parser("verify", help="Cross-check all atoms and subset reconstructions")
    verification.add_argument("input")
    verification.add_argument("--tol", default="1e-8")
    verification.add_argument("--maximum-sources", type=int, choices=range(1, 5), default=3)
    verification.add_argument("--output")
    return command


def main(argv=None):
    command = parser()
    args = command.parse_args(argv)
    try:
        distribution = load_input(args.input)
        if args.command == "query":
            answer = query(
                distribution,
                args.alpha,
                kind=args.kind,
                method=args.method,
                tol=args.tol,
                algorithm=args.algorithm,
                max_seconds=args.max_seconds,
                max_cells=args.max_cells,
                max_level=args.max_level,
                max_terms=args.max_terms,
                maximum_sources=args.maximum_sources,
            )
            completed = answer["status"] == "success" and (
                args.method == "float" or answer["tolerance_met"]
            )
        else:
            from .verification import verify

            answer = verify(distribution, args.tol, maximum_sources=args.maximum_sources)
            completed = answer["status"] == "success"
        payload = json.dumps(serialize(answer), indent=2, allow_nan=False) + "\n"
        if args.output:
            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(payload, encoding="utf-8")
        else:
            print(payload, end="")
        return 0 if completed else 2
    except (OSError, ValueError, TypeError, KeyError, ArithmeticError) as error:
        command.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
