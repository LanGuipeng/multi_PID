from decimal import Decimal
from fractions import Fraction
from time import perf_counter
from .distributions import Distribution
from .lattice import parse_antichain, validate_antichain
from .intervals import Interval, LogExpression


def rational_record(value):
    """Encode exact integers as decimal strings without the JSON digit limit."""
    value = Fraction(value)
    return {
        "numerator": str(Decimal(value.numerator)),
        "denominator": str(Decimal(value.denominator)),
    }


def decode_rational(record):
    fields = [Decimal(record[key]) for key in ("numerator", "denominator")]
    if any(not value.is_finite() or value != value.to_integral_value() for value in fields):
        raise ValueError("Rational fields must be finite integers")
    return Fraction(*(int(value) for value in fields))


def serialize(value):
    if isinstance(value, Fraction):
        return rational_record(value)
    if isinstance(value, Interval):
        return {
            "lower": rational_record(value.lower),
            "upper": rational_record(value.upper),
            "width": rational_record(value.width),
        }
    if isinstance(value, LogExpression):
        return [
            {"coefficient": rational_record(c), "argument": rational_record(x)}
            for c, x in value.terms
        ]
    if isinstance(value, dict):
        return {str(k): serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    return value


def query(P, alpha, kind="atom", method="exact", tol="1e-8", **options):
    """Return a MAR query with its interval, certificate, and computation metrics."""
    start = perf_counter()
    if kind not in ("atom", "cumulative"):
        raise ValueError("kind must be atom or cumulative")
    if isinstance(tol, float):
        raise TypeError("Use a rational tolerance string")
    tol = Fraction(tol)
    if tol <= 0:
        raise ValueError("Positive tolerance required")
    distribution = P if isinstance(P, Distribution) else Distribution(P)
    alpha = (
        parse_antichain(alpha, distribution.n_sources)
        if isinstance(alpha, str) or isinstance(alpha, list)
        else validate_antichain(alpha, distribution.n_sources)
    )
    if method == "exact" or method == "float":
        from .exact import evaluate

        if method == "float":
            options["certify"] = False
    elif method == "endpoint":
        from .endpoint import evaluate
    elif method == "reference":
        if distribution.n_sources > options.get("maximum_sources", 3):
            return {
                "status": "unsupported",
                "reason": "Reference full lattice exceeds configured maximum_sources",
                "method": method,
                "kind": kind,
                "alpha_masks": alpha,
                "interval": None,
                "expression": None,
                "tolerance_met": False,
                "guarantee": "none",
                "timings": {},
                "metrics": {},
            }
        from .reference import evaluate
    else:
        raise ValueError("Unsupported method")
    try:
        result = evaluate(distribution, alpha, kind=kind, tol=tol, **options)
    except TimeoutError as error:
        result = {
            "status": "timeout",
            "reason": str(error),
            "interval": None,
            "expression": None,
            "tolerance_met": False,
            "exact_zero": None,
            "strictly_positive": None,
            "timings": {},
            "metrics": {},
        }
    result.update(
        {"method": method, "kind": kind, "alpha_masks": alpha, "tolerance": tol, "unit": "bit"}
    )
    interval = result.get("interval")
    result["guarantee"] = (
        "deterministic_given_converted_probability_table" if interval is not None else "none"
    )
    if interval is not None:
        result["approximate"] = float((interval.lower + interval.upper) / 2)
        result["width"] = interval.width
    result["structural_status"] = (
        "exact_zero"
        if result.get("exact_zero")
        else ("strictly_positive" if result.get("strictly_positive") else "undetermined")
    )
    result.setdefault("timings", {})["api_total_seconds"] = perf_counter() - start
    result.setdefault("metrics", {}).update(
        {
            "target_alphabet_size": len(distribution.alphabets[0]),
            "source_alphabet_sizes": [len(a) for a in distribution.alphabets[1:]],
            "maximum_denominator_bits": max(
                p.denominator.bit_length() for p in distribution.table.values()
            ),
        }
    )
    result["input_conversion"] = {
        "table_kind": distribution.input_kind,
        "float_policy": distribution.float_policy,
        "normalization": "counts_total" if distribution.input_kind == "counts" else "none",
        "certified_object": "converted_probability_table",
    }
    return result
