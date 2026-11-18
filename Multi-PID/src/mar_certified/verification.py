from fractions import Fraction
from .intervals import LogExpression, certify_expression
from .lattice import antichains, leq
from .exact import build_expression
from .reference import full_lattice


def mutual_information_expression(P, mask):
    coordinates = P.coordinates(mask)
    joint = P.marginal((0,) + coordinates)
    target = P.target_probabilities()
    source = P.marginal(coordinates)
    return LogExpression.from_terms(
        (mass, mass / (target[state[0]] * source[state[1:]]))
        for state, mass in joint.items()
        if mass
    )


def residual_bound(expression, tol):
    interval, terms = certify_expression(expression, Fraction(tol) / 10)
    return max(abs(interval.lower), abs(interval.upper)), interval


def verify(P, tol="1e-8", maximum_sources=3):
    tol = Fraction(tol)
    nodes = antichains(P.n_sources, maximum_sources=maximum_sources)
    values = {a: build_expression(P, a)["expression"] for a in nodes}
    records = []
    for mask in range(1, 1 << P.n_sources):
        total = LogExpression.from_terms(())
        for a, expr in values.items():
            if leq(a, (mask,)):
                total += expr
        difference = total + mutual_information_expression(P, mask).scaled(-1)
        bound, interval = residual_bound(difference, tol)
        records.append(
            {
                "subset": mask,
                "absolute_residual_upper": bound,
                "residual_interval": interval,
                "within_tolerance": bound <= tol,
            }
        )
    independent = []
    reference = full_lattice(P, maximum_sources=maximum_sources)
    for a in nodes:
        bound, interval = residual_bound(values[a] + reference["atoms"][a].scaled(-1), tol)
        independent.append(
            {"alpha": a, "absolute_residual_upper": bound, "within_tolerance": bound <= tol}
        )
    return {
        "status": "success"
        if all(x["within_tolerance"] for x in records + independent)
        else "failed",
        "reconstruction": records,
        "independent_geometry": independent,
        "atoms": len(nodes),
        "interpretation": "Certified absolute residual bounds; interval overlap alone is not a proof",
    }
