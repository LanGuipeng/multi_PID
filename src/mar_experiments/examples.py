"""Taskwise switching and dependent ternary-target decompositions."""

import platform
import sys
import time
from fractions import Fraction as F
from itertools import product
import mpmath as mp
from mar_certified.distributions import Distribution
from mar_certified.exact import evaluate, integral_expression
from mar_certified.reference import full_lattice
from mar_certified.lattice import boundary_sets, leq
from mar_certified.channels import mutual_information_expression
from mar_certified.intervals import LogExpression
from .common import write_json, decimal_output


def switching():
    counts = (8, 12, 32, 48, 76, 4, 19, 1)
    return Distribution(dict(zip(product(range(2), repeat=3), (F(c, 200) for c in counts))))


def ternary(independent=False):
    """T prior (1/2,1/3,1/6); Xi=[T=i-1] XOR Ni.
    (N1,N2) law (14,2,1,3)/20; independent N3~Bernoulli(1/4).
    """
    prior = (F(1, 2), F(1, 3), F(1, 6))
    noise = dict(
        zip(
            product(range(2), repeat=2),
            (F(12, 20), F(4, 20), F(3, 20), F(1, 20))
            if independent
            else (F(14, 20), F(2, 20), F(1, 20), F(3, 20)),
        )
    )
    table = {}
    for t, (n1, n2), n3 in product(range(3), noise, range(2)):
        x = (int(t == 0) ^ n1, int(t == 1) ^ n2, int(t == 2) ^ n3)
        table[(t,) + x] = prior[t] * noise[n1, n2] * (F(3, 4) if n3 == 0 else F(1, 4))
    return Distribution(table)


def mpf(q):
    q = F(q)
    return mp.mpf(q.numerator) / q.denominator


def value(e):
    return mp.fsum((mpf(c) * mp.log(mpf(x), 2) for (c, x) in e.terms))


def factors(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def canonical(e):
    """Prime log basis proves equality of rational log expressions exactly."""
    out = {}
    for c, x in e.terms:
        for sign, n in ((1, x.numerator), (-1, x.denominator)):
            for p, k in factors(n).items():
                out[p] = out.get(p, F()) + sign * k * c
    return {p: c for (p, c) in out.items() if c}


def quadrature(segments):
    result = mp.mpf(0)
    for part in segments:
        p = mpf(part["prior"])
        q = 1 - p
        for left, right, c, d in part["segments"]:
            a = mpf(left)
            b = mp.inf if right is None else mpf(right)
            c = mpf(c)
            d = mpf(d)
            result += p * mp.quad(
                lambda x: (c + d * x) * q * q / (mp.log(2) * x * (p * x + q) ** 2), [a, b]
            )
    return result


def stringify(o):
    if isinstance(o, F):
        return str(o)
    if isinstance(o, mp.mpf):
        return mp.nstr(o, 82)
    if isinstance(o, LogExpression):
        return [[str(c), str(x)] for (c, x) in o.terms]
    if hasattr(o, "lower") and hasattr(o, "upper") and hasattr(o, "width"):
        return {"lower": str(o.lower), "upper": str(o.upper), "width": str(o.width)}
    if isinstance(o, dict):
        return {str(k): stringify(v) for (k, v) in o.items()}
    if isinstance(o, (tuple, list)):
        return [stringify(v) for v in o]
    return o


def check_distribution(d):
    assert sum(d.table.values(), F()) == 1 and min(d.table.values()) >= 0
    for mask in range(1, 1 << d.n_sources):
        coords = d.coordinates(mask)
        for t, p in d.target_probabilities().items():
            if not 0 < p < 1:
                continue
            ch = d.channel(mask, t)
            assert sum((a for (a, b) in ch)) == sum((b for (a, b) in ch)) == 1
            source = d.marginal(coords)
            assert list(source.values()) == [p * a + (1 - p) * b for (a, b) in ch]
    return True


@mp.workdps(80)
def analyze(d, name):
    start = time.perf_counter()
    check_distribution(d)
    d.clear_cache()
    ref = full_lattice(d, maximum_sources=3, max_seconds=90)
    rows = []
    direct = {}
    max_quad_error = mp.mpf(0)
    max_ref_error = mp.mpf(0)
    for alpha in ref["order"]:
        d.clear_cache()
        result = evaluate(d, alpha, kind="atom", tol="1e-15", max_seconds=90)
        expr = result["expression"]
        direct[alpha] = expr
        assert not canonical(expr - ref["atoms"][alpha]), (
            name,
            alpha,
            "symbolic inversion mismatch",
        )
        v = value(expr)
        qv = quadrature(result["segments"])
        assert mpf(result["interval"].lower) <= v <= mpf(result["interval"].upper)
        assert result["interval"].width <= F("1e-15")
        assert result["exact_zero"] or result["interval"].lower > 0
        assert v >= -mp.mpf("1e-75")
        qerr = abs(v - qv)
        rerr = abs(v - value(ref["atoms"][alpha]))
        assert qerr < mp.mpf("1e-65")
        max_quad_error = max(max_quad_error, qerr)
        max_ref_error = max(max_ref_error, rerr)
        target_contributions = {}
        for part in result["segments"]:
            (component, _) = integral_expression(part["prior"], part["segments"])
            target_contributions[part["target"]] = value(component.scaled(part["prior"]))
        rows.append(
            {
                "alpha": alpha,
                "boundary": boundary_sets(alpha, d.n_sources),
                "value_80d": v,
                "target_contributions_80d": target_contributions,
                "structural_zero": result["exact_zero"],
                "interval": result["interval"],
                "width": float(result["interval"].width),
                "log_series_terms": result["series_terms"],
                "expression": expr,
                "prime_log_expression": canonical(expr),
                "segments": result["segments"],
                "metrics": result["metrics"],
                "timings": result["timings"],
                "quadrature_error": qerr,
                "reference_error": rerr,
                "exact_inversion_equality": True,
            }
        )
    mis = []
    for mask in range(1, 1 << d.n_sources):
        expr = mutual_information_expression(d, mask)
        reconstruction = sum(
            (e for (a, e) in direct.items() if leq(a, (mask,))), LogExpression.from_terms(())
        )
        assert not canonical(expr - reconstruction), (name, mask, "reconstruction mismatch")
        coords = d.coordinates(mask)
        joint = d.marginal((0,) + coords)
        marg = d.marginal(coords)
        prior = d.target_probabilities()
        independently = mp.fsum(
            (
                mpf(m) * mp.log(mpf(m) / mpf(prior[s[0]] * marg[s[1:]]), 2)
                for (s, m) in joint.items()
            )
        )
        assert abs(independently - value(expr)) < mp.mpf("1e-70")
        mis.append({"mask": mask, "value_80d": value(expr), "exact_reconstruction_equality": True})
    return {
        "name": name,
        "table": [
            {
                "state": k,
                "probability": v,
                "count_over_480" if name.startswith("ternary") else "count_over_200": int(
                    v * (480 if name.startswith("ternary") else 200)
                ),
            }
            for (k, v) in sorted(d.table.items())
        ],
        "priors": d.target_probabilities(),
        "nodes": len(rows),
        "atoms": rows,
        "mutual_informations": mis,
        "normalization_and_marginal_channels_checked": True,
        "max_quadrature_error": max_quad_error,
        "max_reference_error": max_ref_error,
        "full_geometry_timings": ref["timings"],
        "elapsed_seconds": time.perf_counter() - start,
    }


@mp.workdps(80)
@decimal_output
def run(output_dir):
    report = {
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "mpmath": mp.__version__,
            "mpmath_decimal_digits": 80,
            "input_geometry_arithmetic": "exact Fraction",
            "certification": "rational atanh logarithm series, tolerance 1e-15",
            "clock": "time.perf_counter",
            "per_query_budget_seconds": 90,
        },
        "switching": analyze(switching(), "switching"),
        "ternary": analyze(ternary(), "ternary"),
        "ternary_independent_control": analyze(ternary(True), "ternary_independent_control"),
    }
    s = report["switching"]
    r = next((x for x in s["atoms"] if x["alpha"] == (1, 2)))
    s["redundancy"] = r["value_80d"]
    s["minimum_singleton_mi"] = min(
        (x["value_80d"] for x in s["mutual_informations"] if x["mask"] in (1, 2))
    )
    s["gap"] = s["minimum_singleton_mi"] - s["redundancy"]
    write_json(output_dir, "examples_results.json", stringify(report))
    return report
