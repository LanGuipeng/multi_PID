import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import pytest

from mar_certified import Distribution, query
from mar_certified.api import decode_rational, serialize
from mar_certified.intervals import LogExpression
from mar_certified.lattice import boundary_sets
from mar_certified.verification import verify


@pytest.fixture
def switching():
    counts = [8, 12, 32, 48, 76, 4, 19, 1]
    return Distribution(dict(zip(product((0, 1), repeat=3), counts)), kind="counts")


def test_copy_and_xor_have_distinct_atomic_support():
    copy = Distribution({(0, 0, 0): 1, (1, 1, 1): 1}, kind="counts")
    xor = Distribution({(a ^ b, a, b): 1 for a, b in product((0, 1), repeat=2)}, kind="counts")
    for distribution, nonzero in [(copy, "1|2"), (xor, "1,2")]:
        for alpha in ["1|2", "1", "2", "1,2"]:
            result = query(distribution, alpha, tol="1e-20")
            assert result["tolerance_met"]
            assert result["interval"].contains(int(alpha == nonzero))
            if alpha != nonzero:
                assert result["exact_zero"]


def test_switching_matches_analytic_expression(switching):
    expected = LogExpression.from_terms(
        [
            (Fraction(11, 5), 2),
            (Fraction(-1, 2), 5),
            (Fraction(-21, 10), 3),
            (Fraction(7, 40), 7),
            (Fraction(19, 40), 19),
        ]
    ).interval(32)
    result = query(switching, [[1], [2]], tol="1e-12")
    assert result["interval"].lower <= expected.lower <= expected.upper <= result["interval"].upper
    assert result["approximate"] == pytest.approx(0.2196673812876747, abs=1e-13)
    for algorithm in ["enumerate", "optimized"]:
        check = query(switching, "1|2", algorithm=algorithm, tol="1e-12")
        assert check["expression"] == result["expression"]
    minimum_mi = min(query(switching, a, kind="cumulative")["approximate"] for a in ["1", "2"])
    assert minimum_mi - result["approximate"] == pytest.approx(0.05840452382496298, abs=1e-8)


def test_nonuniform_dependent_target_reconstructs_all_subsets():
    from mar_experiments.examples import ternary

    result = verify(ternary(), tol="1e-15", maximum_sources=3)
    assert result["status"] == "success"
    assert result["atoms"] == 18
    assert len(result["reconstruction"]) == 7
    assert all(row["within_tolerance"] for row in result["independent_geometry"])


def test_small_positive_atom_retains_a_positive_certificate():
    delta = Fraction(1, 10**12)
    table = {
        (2 * u + v, 2 * u + w, u): (1 + (2 * delta if w == v else -2 * delta)) / 8
        for u, v, w in product((0, 1), repeat=3)
    }
    result = query(Distribution(table), "1", tol="1e-32")
    assert result["interval"].lower > 0
    assert result["interval"].width <= Fraction("1e-32")
    assert result["approximate"] == pytest.approx(2.885390081777927e-24, rel=1e-8, abs=0)


def test_endpoint_budget_preserves_valid_interval(switching):
    exact = query(switching, "1|2", tol="1e-30")["interval"]
    result = query(switching, "1|2", method="endpoint", tol="1e-12", max_cells=0)
    assert not result["tolerance_met"]
    assert result["status"] != "success"
    assert result["interval"].lower <= exact.lower <= exact.upper <= result["interval"].upper
    reached = query(switching, "1|2", method="endpoint", tol="1e-1", max_cells=2048, max_seconds=10)
    assert reached["tolerance_met"]
    assert reached["interval"].width <= Fraction("1e-1")
    assert reached["interval"].lower <= exact.lower <= exact.upper <= reached["interval"].upper


def test_probability_and_antichain_validation():
    with pytest.raises(TypeError):
        Distribution({(0, 0): 0.5, (1, 1): 0.5})
    with pytest.raises(ValueError):
        Distribution({(0, 0): "1/3"})
    assert boundary_sets((1, 6), 3) == (2, 4)
    assert boundary_sets((1, 6), 4) == (10, 12)
    with pytest.raises(ValueError):
        query({(0, 0, 0): 1}, [[1], [1, 2]])


def test_zero_prior_constant_target_and_singular_columns():
    copy = Distribution({(0, 0, 0): "1/3", (1, 1, 1): "2/3", (2, 0, 0): 0})
    result = query(copy, "1|2", tol="1e-15")
    assert result["approximate"] == pytest.approx(0.9182958340544896, abs=1e-14)
    assert any(part["segments"][-1][1] is None for part in result["segments"])
    constant = query({(0, 0, 0): Fraction(1, 2), (0, 1, 1): Fraction(1, 2)}, "1")
    assert constant["interval"].lower == constant["interval"].upper == 0
    assert constant["exact_zero"]


def test_certificates_survive_json_round_trip(switching):
    result = query(switching, "1|2", tol="1e-20")
    decoded = json.loads(json.dumps(serialize(result)))
    lower = decode_rational(decoded["interval"]["lower"])
    upper = decode_rational(decoded["interval"]["upper"])
    assert lower == result["interval"].lower
    assert upper == result["interval"].upper
    assert upper - lower <= Fraction("1e-20")


def test_reference_records_preserve_budget_failures():
    path = Path(__file__).resolve().parents[1] / "results/reference/precision_results.json"
    rows = json.loads(path.read_text(encoding="utf-8"))["runs"]
    for row in rows:
        assert Fraction(row["exact"]["width"]) <= Fraction(row["epsilon"])
        if "endpoint" in row:
            endpoint = row["endpoint"]
            assert endpoint["contains_exact_certificate"]
            assert endpoint["tolerance_met"] == (
                Fraction(endpoint["interval"]["width"]) <= Fraction(row["epsilon"])
            )


def test_large_rational_serialization_without_global_digit_changes():
    import sys

    value = Fraction(10**5000 + 1, 10**5000 + 3)
    before = sys.get_int_max_str_digits() if hasattr(sys, "get_int_max_str_digits") else None
    encoded = json.loads(json.dumps(serialize(value)))
    assert decode_rational(encoded) == value
    if before is not None:
        assert sys.get_int_max_str_digits() == before


@pytest.mark.parametrize("invalid", ["1.9", "NaN", "Infinity", "-Infinity"])
def test_rational_decoder_rejects_noninteger_fields(invalid):
    with pytest.raises(ValueError, match="finite integers"):
        decode_rational({"numerator": invalid, "denominator": "2"})
