from fractions import Fraction
from hashlib import sha256
import json


def rational(value, float_policy=None):
    if isinstance(value, bool):
        raise TypeError("Boolean is not a probability")
    if isinstance(value, float):
        if float_policy == "decimal_string":
            return Fraction(str(value))
        if float_policy == "binary_exact":
            return Fraction.from_float(value)
        raise TypeError("Native float requires decimal_string or binary_exact policy")
    if not isinstance(value, (int, str, Fraction)):
        raise TypeError("Use integer, Fraction or rational string")
    return Fraction(value)


class Distribution:
    """Exactly normalized joint probabilities indexed by (target, source_1, ...)."""

    def __init__(self, table, kind="probabilities", float_policy=None, cache=True):
        if kind not in ("probabilities", "counts"):
            raise ValueError("Unknown table kind")
        if float_policy not in (None, "decimal_string", "binary_exact"):
            raise ValueError("Unknown float conversion policy")
        if not table:
            raise ValueError("Empty table")
        if any(not isinstance(state, tuple) for state in table):
            raise TypeError("State keys must be tuples")
        lengths = {len(state) for state in table}
        if len(lengths) != 1 or next(iter(lengths)) < 2:
            raise ValueError("States must be (T,X1,...,Xn), n >= 1")
        for state in table:
            for value in state:
                hash(value)
        converted = {state: rational(mass, float_policy) for state, mass in table.items()}
        if any(mass < 0 for mass in converted.values()):
            raise ValueError("Negative mass")
        total = sum(converted.values(), Fraction())
        if kind == "counts":
            if any(mass.denominator != 1 for mass in converted.values()) or total <= 0:
                raise ValueError("Counts must be nonnegative integers with positive total")
            converted = {state: mass / total for state, mass in converted.items()}
        elif total != 1:
            raise ValueError(
                "Probability masses must sum exactly to one; no implicit normalization"
            )
        self.n_sources = next(iter(lengths)) - 1
        self.table = {state: mass for state, mass in converted.items() if mass}
        self.alphabets = tuple(
            tuple(dict.fromkeys(state[i] for state in table)) for i in range(self.n_sources + 1)
        )
        self.input_kind = kind
        self.float_policy = float_policy
        self.cache_enabled = bool(cache)
        self._marginals = {}
        self.cache_hits = 0
        self.cache_misses = 0

    @classmethod
    def from_rows(cls, rows, kind="probabilities", float_policy=None, cache=True):
        table = {}
        for row in rows:
            state = tuple(row["state"])
            if state in table:
                raise ValueError("Duplicate state row")
            table[state] = row.get("mass", row.get("count"))
        return cls(table, kind, float_policy, cache)

    def marginal(self, coordinates):
        coordinates = tuple(coordinates)
        if not coordinates or len(set(coordinates)) != len(coordinates):
            raise ValueError("Coordinates must be nonempty and distinct")
        if any(type(i) is not int or not 0 <= i <= self.n_sources for i in coordinates):
            raise ValueError("Invalid coordinate")
        if self.cache_enabled and coordinates in self._marginals:
            self.cache_hits += 1
            return dict(self._marginals[coordinates])
        self.cache_misses += 1
        values = {}
        for state, mass in self.table.items():
            key = tuple(state[i] for i in coordinates)
            values[key] = values.get(key, Fraction()) + mass
        if self.cache_enabled:
            self._marginals[coordinates] = values
        return dict(values)

    def target_probabilities(self):
        observed = {key[0]: mass for key, mass in self.marginal((0,)).items()}
        return {target: observed.get(target, Fraction()) for target in self.alphabets[0]}

    def coordinates(self, mask):
        if type(mask) is not int or not 0 < mask < 1 << self.n_sources:
            raise ValueError("Invalid source mask")
        return tuple(i + 1 for i in range(self.n_sources) if mask & (1 << i))

    def channel(self, mask, target):
        indices = self.coordinates(mask)
        prior = self.target_probabilities().get(target, Fraction())
        if not 0 < prior < 1:
            raise ValueError("Conditional channel requires 0 < prior < 1")
        source = self.marginal(indices)
        joint = self.marginal((0,) + indices)
        return tuple(
            (
                joint.get((target,) + state, Fraction()) / prior,
                (mass - joint.get((target,) + state, Fraction())) / (1 - prior),
            )
            for state, mass in source.items()
        )

    def fingerprint(self):
        rows = sorted(
            ((repr(state), mass.numerator, mass.denominator) for state, mass in self.table.items())
        )
        return sha256(json.dumps(rows, ensure_ascii=False).encode()).hexdigest()

    def clear_cache(self):
        self._marginals.clear()
        self.cache_hits = 0
        self.cache_misses = 0
