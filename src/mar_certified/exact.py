from fractions import Fraction
from math import log2
from time import perf_counter
from .intervals import Interval, LogExpression, certify_expression
from .lattice import validate_antichain, boundary_sets
from .piecewise import decision_segments


def primitive_terms(prior, intercept, slope, point):
    q = 1 - prior
    if point is None:
        return ((intercept, 1 / prior),), q * slope / prior
    if point == 0:
        if intercept:
            raise ArithmeticError("Nonintegrable primitive at zero")
        return (), Fraction()
    denominator = prior * point + q
    return ((intercept, point / denominator),), q * (intercept + slope * point) / denominator


def integral_expression(prior, segments):
    prior = Fraction(prior)
    if not 0 < prior < 1:
        raise ValueError("Nondegenerate prior required")
    terms = []
    rational_residual = Fraction()
    previous_right = Fraction()
    previous_value = Fraction()
    endpoint_terms = []
    tail_seen = False
    for left, right, intercept, slope in segments:
        if tail_seen:
            raise ArithmeticError("A tail segment must be last")
        tail_seen = right is None
        if left < previous_right or (right is not None and right <= left):
            raise ArithmeticError("Segments must be ordered and nonempty")
        start_value = intercept + slope * left
        if left == previous_right:
            if start_value != previous_value:
                raise ArithmeticError("Discontinuous complete envelope")
        elif previous_value != 0 or start_value != 0:
            raise ArithmeticError("Gap endpoints must have zero envelope values")
        if right is None and slope:
            raise ArithmeticError("Actual finite-channel tail must have zero slope")
        lower_terms, lower_residual = primitive_terms(prior, intercept, slope, left)
        upper_terms, upper_residual = primitive_terms(prior, intercept, slope, right)
        terms.extend((-c, x) for c, x in lower_terms)
        terms.extend(upper_terms)
        rational_residual += upper_residual - lower_residual
        endpoint_terms.append((lower_residual, upper_residual))
        if right is not None:
            previous_right = right
            previous_value = intercept + slope * right
        else:
            previous_value = Fraction()
    if previous_value != 0:
        raise ArithmeticError("Finite positive envelope cannot end at a positive value")
    if rational_residual:
        raise ArithmeticError("Primitive rational endpoint terms did not cancel")
    return (
        LogExpression.from_terms(terms),
        {"rational_residual": rational_residual, "endpoint_terms": endpoint_terms},
    )


def build_expression(
    distribution,
    alpha,
    kind="atom",
    algorithm="enumerate",
    max_seconds=60,
    max_segments=100000,
    max_boundary_subsets=1048576,
    **options,
):
    if kind not in ("atom", "cumulative"):
        raise ValueError("Unknown query kind")
    started = perf_counter()
    deadline = started + max_seconds if max_seconds is not None else float("inf")

    def check():
        if perf_counter() > deadline:
            raise TimeoutError("Exact query time budget exceeded")

    alpha = validate_antichain(alpha, distribution.n_sources)
    priors = distribution.target_probabilities()
    timings = {
        "boundary": 0.0,
        "marginalization": 0.0,
        "envelope": 0.0,
        "integration": 0.0,
        "certification": 0.0,
    }
    expression = LogExpression.from_terms(())
    metrics = {
        "n_sources": distribution.n_sources,
        "support_size": len(distribution.table),
        "alpha_size": len(alpha),
        "boundary_size": 0,
        "channel_columns": 0,
        "base_intervals": 0,
        "effective_segments": 0,
        "probability_bit_length": sum(
            p.numerator.bit_length() + p.denominator.bit_length()
            for p in distribution.table.values()
        ),
        "state_encoding_bytes": sum(len(repr(state).encode()) for state in distribution.table),
        "maximum_denominator_bits": max(
            p.denominator.bit_length() for p in distribution.table.values()
        ),
        "target_alphabet_size": len(distribution.alphabets[0]),
        "source_alphabet_sizes": tuple(len(a) for a in distribution.alphabets[1:]),
    }
    metrics["input_bit_length"] = (
        metrics["probability_bit_length"] + 8 * metrics["state_encoding_bytes"]
    )
    if sum(p > 0 for p in priors.values()) <= 1:
        return {
            "expression": expression,
            "segments": [],
            "endpoint_audit": [],
            "exact_zero": True,
            "strictly_positive": False,
            "timings": timings,
            "metrics": metrics,
            "status": "success",
            "deadline": deadline,
        }
    stage = perf_counter()
    competitors = (
        boundary_sets(alpha, distribution.n_sources, max_boundary_subsets, check)
        if kind == "atom"
        else ()
    )
    timings["boundary"] = perf_counter() - stage
    metrics["boundary_size"] = len(competitors)
    metrics["query_masks"] = alpha
    metrics["boundary_masks"] = competitors
    all_segments = []
    audits = []
    for target, prior in priors.items():
        if not 0 < prior < 1:
            continue
        check()
        stage = perf_counter()
        channels = tuple(distribution.channel(mask, target) for mask in alpha + competitors)
        timings["marginalization"] += perf_counter() - stage
        metrics["channel_columns"] += sum(len(channel) for channel in channels)
        metrics["base_intervals"] += len(
            {Fraction(), Fraction(1)}
            | {a / b for channel in channels for a, b in channel if a and b}
        )
        stage = perf_counter()
        segments = decision_segments(
            channels, len(alpha), kind == "atom", algorithm, check, max_segments
        )
        timings["envelope"] += perf_counter() - stage
        metrics["effective_segments"] += len(segments)
        if metrics["effective_segments"] > max_segments:
            raise TimeoutError("Total segment budget exceeded")
        stage = perf_counter()
        component, audit = integral_expression(prior, segments)
        expression += component.scaled(prior)
        timings["integration"] += perf_counter() - stage
        all_segments.append({"target": target, "prior": prior, "segments": segments})
        audits.append(audit)
    metrics["cache_hits"] = distribution.cache_hits
    metrics["cache_misses"] = distribution.cache_misses
    metrics["expression_terms"] = len(expression.terms)
    zero = metrics["effective_segments"] == 0
    return {
        "expression": expression,
        "segments": all_segments,
        "endpoint_audit": audits,
        "exact_zero": zero,
        "strictly_positive": not zero,
        "timings": timings,
        "metrics": metrics,
        "status": "success",
        "deadline": deadline,
    }


def evaluate(
    distribution, alpha, kind="atom", tol="1e-8", algorithm="enumerate", certify=True, **budgets
):
    started = perf_counter()
    epsilon = Fraction(tol)
    if epsilon <= 0:
        raise ValueError("Positive tolerance required")
    result = build_expression(distribution, alpha, kind, algorithm, **budgets)
    deadline = result.pop("deadline")
    if not certify:
        result.update(
            {
                "interval": None,
                "approximate": sum(
                    float(c) * log2(float(x)) for c, x in result["expression"].terms
                ),
                "tolerance_met": False,
                "guarantee": "floating_evaluation_of_exact_expression",
                "series_terms": 0,
            }
        )
    else:
        stage = perf_counter()
        if result["exact_zero"]:
            interval, terms = Interval(0, 0), 0
        else:
            interval, terms = certify_expression(
                result["expression"], epsilon, max_terms=budgets.get("max_terms"), deadline=deadline
            )
        result["timings"]["certification"] = perf_counter() - stage
        result.update(
            {
                "interval": interval,
                "approximate": float((interval.lower + interval.upper) / 2),
                "tolerance_met": interval.width <= epsilon,
                "guarantee": "deterministic_rational_enclosure",
                "series_terms": terms,
            }
        )
    result["timings"]["total"] = perf_counter() - started
    return result
