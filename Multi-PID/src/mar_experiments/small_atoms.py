"""Certified small atoms, floating-point cancellation, and support cases."""

import math
import platform
import sys
from fractions import Fraction as F
import mpmath as mp
from mar_certified.distributions import Distribution
from mar_certified.exact import evaluate, build_expression
from mar_certified.endpoint import evaluate as endpoint
from mar_certified.reference import full_lattice
from mar_certified.lattice import leq, antichains
from mar_certified.intervals import LogExpression, Interval, clear_caches, certify_expression
from mar_certified.verification import verify
from .common import write_json, decimal_output


def mpr(x):
    x = F(x)
    return mp.mpf(x.numerator) / x.denominator


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, Interval):
        return {
            "lower": str(x.lower),
            "upper": str(x.upper),
            "lower_float": float(x.lower),
            "upper_float": float(x.upper),
            "width_float": float(x.width),
        }
    if isinstance(x, LogExpression):
        return {"terms": [[str(c), str(a)] for (c, a) in x.terms]}
    if isinstance(x, dict):
        return {str(k): encode(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def small_distribution(delta):
    delta = F(delta)
    assert 0 <= delta < F(1, 2)
    table = {}
    for u in (0, 1):
        for v in (0, 1):
            for noise in (0, 1):
                table[(u, v), (u, v ^ noise), u] = (
                    F(1, 2) + delta if noise == 0 else F(1, 2) - delta
                ) / 4
    return Distribution(table)


def analytic_expression(delta):
    delta = F(delta)
    return LogExpression.from_terms(
        ((F(1, 2) + delta, 1 + 2 * delta), (F(1, 2) - delta, 1 - 2 * delta))
    )


def float_mobius(P):
    full = full_lattice(P)
    atoms = {}
    for alpha in full["order"]:
        value = sum((float(c) * math.log2(float(x)) for (c, x) in full["cumulative"][alpha].terms))
        atoms[alpha] = value - sum((v for (b, v) in atoms.items() if leq(b, alpha)))
    return atoms[1,]


def validate_channels(P):
    assert sum(P.table.values()) == 1
    checks = 0
    for target, prior in P.target_probabilities().items():
        if 0 < prior < 1:
            for mask in range(1, 1 << P.n_sources):
                cols = P.channel(mask, target)
                assert sum((a for (a, b) in cols)) == 1 and sum((b for (a, b) in cols)) == 1
                checks += 1
    return checks


def validate_nonnegative(P):
    """Certify nonnegative affine envelopes, including the constant tail."""
    report = []
    for alpha in antichains(P.n_sources):
        result = build_expression(P, alpha)
        for target in result["segments"]:
            for left, right, intercept, slope in target["segments"]:
                assert intercept + slope * left >= 0
                if right is None:
                    assert slope == 0
                else:
                    assert intercept + slope * right >= 0
        report.append(
            {"alpha": alpha, "structural_zero": result["exact_zero"], "nonnegative_segments": True}
        )
    return report


def compact_exact(result):
    return {k: v for (k, v) in result.items() if k not in ("endpoint_audit",)}


@mp.workdps(100)
@decimal_output
def run(output_dir):
    output = {
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "mpmath": mp.__version__,
            "reference_decimal_digits": mp.mp.dps,
            "timer": "time.perf_counter",
            "timing_policy": "single serial run; logarithm and marginal caches cleared before exact query",
        },
        "family": [],
        "endpoint": [],
        "support": [],
    }
    for k in (1, 3, 6, 9, 12):
        delta = F(1, 10**k)
        P = small_distribution(delta)
        clear_caches()
        eps = delta**2 / F(10**8)
        exact = evaluate(P, (1,), tol=eps, max_seconds=30)
        expr = analytic_expression(delta)
        ref = sum((mpr(c) * mp.log(mpr(a), 2) for (c, a) in expr.terms))
        interval = exact["interval"]
        assert mpr(interval.lower) <= ref <= mpr(interval.upper)
        assert interval.width <= eps and interval.lower > 0
        difference = exact["expression"] - expr
        (residual, _) = certify_expression(difference, eps / F(100))
        assert max(abs(residual.lower), abs(residual.upper)) <= eps / F(100)
        naive = float_mobius(P)
        direct_float = sum((float(c) * math.log2(float(a)) for (c, a) in exact["expression"].terms))
        checked = verify(P, tol=eps / F(100))
        assert checked["status"] == "success"
        record = {
            "k": k,
            "delta": delta,
            "epsilon": eps,
            "reference": mp.nstr(ref, 70),
            "exact_midpoint_absolute_error": mp.nstr(
                abs(mpr((interval.lower + interval.upper) / 2) - ref), 30
            ),
            "float_mobius": naive,
            "float_mobius_absolute_error": mp.nstr(abs(mp.mpf(naive) - ref), 30),
            "float_mobius_relative_error": mp.nstr(abs(mp.mpf(naive) - ref) / ref, 30),
            "float_direct_expression": direct_float,
            "float_direct_absolute_error": mp.nstr(abs(mp.mpf(direct_float) - ref), 30),
            "channel_normalizations": validate_channels(P),
            "verification": checked,
            "all_atom_nonnegativity": validate_nonnegative(P),
            "analytic_residual": residual,
            "exact": compact_exact(exact),
            "lambda_minus": 3 * (1 - 2 * delta) / (1 + 2 * delta),
            "lambda_plus": 3 * (1 + 2 * delta) / (1 - 2 * delta),
            "peak": 2 * delta,
            "reference_inside_interval": True,
        }
        output["family"].append(record)
    for delta, tol in (
        (F(1, 10), F(1, 10)),
        (F(1, 1000), F(1, 100)),
        (F(1, 1000), F(1, 10**8)),
        (F(0), F(1, 100)),
    ):
        clear_caches()
        P = small_distribution(delta)
        result = endpoint(
            P, (1,), tol=tol, max_seconds=12, max_cells=16384, max_level=12, max_terms=128
        )
        expr = analytic_expression(delta)
        ref = sum((mpr(c) * mp.log(mpr(a), 2) for (c, a) in expr.terms))
        assert mpr(result["interval"].lower) <= ref <= mpr(result["interval"].upper)
        (independent_certificate, _) = certify_expression(expr, F(1, 10**40))
        assert result["interval"].lower <= independent_certificate.lower
        assert independent_certificate.upper <= result["interval"].upper
        output["endpoint"].append(
            {
                "delta": delta,
                "epsilon": tol,
                "max_seconds": 12,
                "max_cells": 16384,
                "reference": mp.nstr(ref, 70),
                "reference_inside_interval": True,
                "contains_analytic_rational_certificate": True,
                "analytic_rational_certificate": independent_certificate,
                "result": result,
            }
        )
    cases = [
        (
            "zero_probability_target",
            Distribution({(0, 0, 0): F(1, 2), (1, 1, 1): F(1, 2), (2, 0, 0): 0}),
            (1, 2),
        ),
        ("constant_target", Distribution({(0, 0, 0): F(1, 2), (0, 1, 1): F(1, 2)}), (1,)),
        (
            "singular_observation_column",
            Distribution({(0, 0, 0): F(1, 3), (1, 1, 1): F(2, 3)}),
            (1, 2),
        ),
        (
            "rare_target_prior",
            Distribution({(0, 0, 0): F(1, 10**6), (1, 1, 1): 1 - F(1, 10**6)}),
            (1, 2),
        ),
        (
            "identical_observation_unique_atom",
            Distribution({(0, 0, 0): F(1, 3), (1, 1, 1): F(2, 3)}),
            (1,),
        ),
    ]
    for name, P, alpha in cases:
        clear_caches()
        exact = evaluate(P, alpha, tol=F(1, 10**16), max_seconds=30)
        checked = verify(P, tol=F(1, 10**16))
        assert checked["status"] == "success"
        tail = []
        for target, prior in P.target_probabilities().items():
            if 0 < prior < 1:
                cols = P.channel(alpha[0], target)
                tail.append(
                    {
                        "target": target,
                        "prior": prior,
                        "sum_a_for_b_zero": sum((a for (a, b) in cols if not b), F()),
                    }
                )
        output["support"].append(
            {
                "name": name,
                "table": [{"state": s, "probability": p} for (s, p) in P.table.items()],
                "declared_target_alphabet": P.alphabets[0],
                "alpha": alpha,
                "exact": compact_exact(exact),
                "verification": checked,
                "singular_columns": tail,
                "channel_normalizations": validate_channels(P),
                "all_atom_nonnegativity": validate_nonnegative(P),
            }
        )
    write_json(output_dir, "small_atoms_results.json", encode(output))
    return output
