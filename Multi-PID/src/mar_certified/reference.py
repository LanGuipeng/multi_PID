from fractions import Fraction
from itertools import combinations
from .intervals import LogExpression, certify_expression
from time import perf_counter
from math import log2
from .lattice import leq, validate_antichain, topological_antichains


def cross(origin, first, second):
    return (first[0] - origin[0]) * (second[1] - origin[1]) - (first[1] - origin[1]) * (
        second[0] - origin[0]
    )


def zonogon_vertices(columns):
    ordered = sorted(
        columns, key=lambda column: (column[0] == 0, column[1] / column[0] if column[0] else 0)
    )
    path = [(Fraction(), Fraction())]
    for a, b in ordered:
        x, y = path[-1]
        path.append((x + a, y + b))
    return tuple(path + [(1 - x, 1 - y) for x, y in path])


def join_columns(channels):
    points = sorted({point for channel in channels for point in zonogon_vertices(channel)})
    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    columns = tuple((v[0] - u[0], v[1] - u[1]) for u, v in zip(lower, lower[1:]))
    if any(a < 0 or b < 0 for a, b in columns):
        raise ArithmeticError("Invalid join boundary")
    return columns


def geometry_cumulative_expression(distribution, alpha):
    alpha = validate_antichain(alpha, distribution.n_sources)
    expression = LogExpression.from_terms(())
    for target, prior in distribution.target_probabilities().items():
        if not 0 < prior < 1:
            continue
        channels = tuple(distribution.channel(mask, target) for mask in alpha)
        for size in range(1, len(channels) + 1):
            sign = (-1) ** (size - 1)
            for chosen in combinations(channels, size):
                columns = join_columns(chosen)
                expression += LogExpression.from_terms(
                    (prior * sign * a, a / (prior * a + (1 - prior) * b)) for a, b in columns if a
                )
    return expression


def full_lattice(distribution, maximum_sources=3, max_seconds=60):
    started = perf_counter()
    deadline = started + max_seconds if max_seconds is not None else float("inf")
    order = topological_antichains(distribution.n_sources, maximum_sources)
    cumulative = {}
    atoms = {}
    geometry_time = 0.0
    inversion_time = 0.0
    for alpha in order:
        if perf_counter() > deadline:
            raise TimeoutError("Reference full lattice time budget exceeded")
        stage = perf_counter()
        cumulative[alpha] = geometry_cumulative_expression(distribution, alpha)
        geometry_time += perf_counter() - stage
        stage = perf_counter()
        expression = cumulative[alpha]
        for beta, atom in atoms.items():
            if leq(beta, alpha):
                expression += atom.scaled(-1)
        atoms[alpha] = expression
        inversion_time += perf_counter() - stage
    return {
        "atoms": atoms,
        "cumulative": cumulative,
        "order": order,
        "timings": {
            "geometry": geometry_time,
            "inversion": inversion_time,
            "total": perf_counter() - started,
        },
    }


def reference_atom_expression(distribution, alpha, maximum_sources=3):
    return full_lattice(distribution, maximum_sources)["atoms"][
        validate_antichain(alpha, distribution.n_sources)
    ]


def evaluate(
    distribution,
    alpha,
    kind="atom",
    tol="1e-8",
    maximum_sources=3,
    max_seconds=60,
    certify=True,
    max_terms=None,
    **budgets,
):
    started = perf_counter()
    epsilon = Fraction(tol)
    if epsilon <= 0:
        raise ValueError("Positive tolerance required")
    alpha = validate_antichain(alpha, distribution.n_sources)
    if kind == "atom":
        full = full_lattice(distribution, maximum_sources, max_seconds)
        expression = full["atoms"][alpha]
        timings = full["timings"]
        nodes = len(full["atoms"])
    elif kind == "cumulative":
        if distribution.n_sources > maximum_sources:
            raise ValueError("Reference exceeds configured source budget")
        expression = geometry_cumulative_expression(distribution, alpha)
        timings = {"geometry": perf_counter() - started, "inversion": 0.0}
        nodes = 1
    else:
        raise ValueError("Unknown query kind")
    deadline = started + max_seconds if max_seconds is not None else float("inf")
    stage = perf_counter()
    if certify:
        interval, terms = certify_expression(
            expression, epsilon, max_terms=max_terms, deadline=deadline
        )
        approximate = float((interval.lower + interval.upper) / 2)
    else:
        interval, terms = None, 0
        approximate = sum(float(c) * log2(float(x)) for c, x in expression.terms)
    timings["certification"] = perf_counter() - stage
    timings["total"] = perf_counter() - started
    return {
        "expression": expression,
        "interval": interval,
        "approximate": approximate,
        "tolerance_met": interval is not None and interval.width <= epsilon,
        "exact_zero": not expression.terms,
        "strictly_positive": interval is not None and interval.lower > 0,
        "guarantee": "deterministic_rational_enclosure"
        if certify
        else "floating_evaluation_of_exact_expression",
        "series_terms": terms,
        "timings": timings,
        "status": "success",
        "metrics": {
            "n_sources": distribution.n_sources,
            "support_size": len(distribution.table),
            "lattice_nodes": nodes,
            "expression_terms": len(expression.terms),
        },
    }
