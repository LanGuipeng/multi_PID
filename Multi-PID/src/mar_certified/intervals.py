from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from time import perf_counter


@dataclass(frozen=True)
class Interval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self):
        object.__setattr__(self, "lower", Fraction(self.lower))
        object.__setattr__(self, "upper", Fraction(self.upper))
        if self.lower > self.upper:
            raise ValueError("Reversed interval")

    @classmethod
    def point(cls, value):
        return cls(value, value)

    @property
    def width(self):
        return self.upper - self.lower

    def contains(self, value):
        return self.lower <= value <= self.upper

    def positive_part(self):
        return Interval(max(Fraction(), self.lower), max(Fraction(), self.upper))

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval.point(other)
        return Interval(self.lower + other.lower, self.upper + other.upper)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Interval) else -Fraction(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval.point(other)
        products = tuple(
            a * b for a in (self.lower, self.upper) for b in (other.lower, other.upper)
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lower <= 0 <= self.upper:
            raise ZeroDivisionError("Interval contains zero")
        return Interval(1 / self.upper, 1 / self.lower)

    def __truediv__(self, other):
        other = other if isinstance(other, Interval) else Interval.point(other)
        return self * other.reciprocal()

    def __rtruediv__(self, other):
        return self.reciprocal() * other


def interval_min(values):
    values = tuple(values)
    if not values:
        raise ValueError("Minimum of empty interval family")
    return Interval(min(v.lower for v in values), min(v.upper for v in values))


def interval_max(values):
    values = tuple(values)
    if not values:
        return Interval.point(0)
    return Interval(max(v.lower for v in values), max(v.upper for v in values))


@lru_cache(maxsize=8192)
def natural_log_unit(value, terms):
    value = Fraction(value)
    if not 1 <= value <= 2 or not isinstance(terms, int) or terms < 1:
        raise ValueError("Require 1 <= value <= 2 and positive integer terms")
    z = (value - 1) / (value + 1)
    power = z
    total = Fraction()
    for index in range(terms):
        total += 2 * power / (2 * index + 1)
        power *= z * z
    return Interval(total, total + 2 * power / ((2 * terms + 1) * (1 - z * z)))


def normalized_log_argument(value):
    value = Fraction(value)
    if value <= 0:
        raise ValueError("Logarithm requires a positive rational")
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    normalized = value / Fraction(2) ** exponent
    while normalized < 1:
        normalized *= 2
        exponent -= 1
    while normalized >= 2:
        normalized /= 2
        exponent += 1
    return exponent, normalized


@lru_cache(maxsize=8192)
def log2_interval(value, terms=16):
    if not isinstance(terms, int) or terms < 1:
        raise ValueError("Positive integer terms required")
    exponent, normalized = normalized_log_argument(value)
    if normalized == 1:
        return Interval.point(exponent)
    return exponent + natural_log_unit(normalized, terms) / natural_log_unit(Fraction(2), terms)


@lru_cache(maxsize=8192)
def natural_log_interval(value, terms=16):
    exponent, normalized = normalized_log_argument(value)
    return exponent * natural_log_unit(Fraction(2), terms) + natural_log_unit(normalized, terms)


def phi_prime_interval(prior, value, terms=16):
    prior, value = Fraction(prior), Fraction(value)
    if not 0 < prior < 1 or value <= 0:
        raise ValueError("Require 0 < prior < 1 and value > 0")
    denominator = prior * value + 1 - prior
    return log2_interval(value / denominator, terms) + Interval.point(
        (1 - prior) / denominator
    ) / natural_log_unit(Fraction(2), terms)


def weight_interval(prior, left, right, terms=16):
    prior, left, right = Fraction(prior), Fraction(left), Fraction(right)
    if not 0 < prior < 1 or not 0 < left <= right:
        raise ValueError("Require valid prior and positive ordered endpoints")
    if left == right:
        return Interval.point(0)
    q = 1 - prior
    first, second = prior * left + q, prior * right + q
    ratio = right * first / (left * second)
    rational_term = q * (1 / second - 1 / first)
    result = log2_interval(ratio, terms) + Interval.point(rational_term) / natural_log_unit(
        Fraction(2), terms
    )
    if result.upper < 0:
        raise ArithmeticError("Negative upper bound for a positive weight")
    return Interval(max(Fraction(), result.lower), result.upper)


@dataclass(frozen=True)
class LogExpression:
    terms: tuple

    @classmethod
    def from_terms(cls, terms):
        coefficients = {}
        for coefficient, argument in terms:
            coefficient, argument = Fraction(coefficient), Fraction(argument)
            if argument <= 0:
                raise ValueError("Nonpositive logarithm argument")
            if argument != 1:
                coefficients[argument] = coefficients.get(argument, Fraction()) + coefficient
        return cls(
            tuple(
                (coefficient, argument)
                for argument, coefficient in sorted(coefficients.items())
                if coefficient
            )
        )

    def interval(self, terms):
        return sum(
            (coefficient * log2_interval(argument, terms) for coefficient, argument in self.terms),
            Interval.point(0),
        )

    def scaled(self, coefficient):
        return LogExpression.from_terms(
            (Fraction(coefficient) * value, argument) for value, argument in self.terms
        )

    def __add__(self, other):
        return LogExpression.from_terms(self.terms + other.terms)

    def __neg__(self):
        return self.scaled(-1)

    def __sub__(self, other):
        return self + (-other)


def certify_expression(expression, epsilon, max_terms=None, deadline=None):
    epsilon = Fraction(epsilon)
    if epsilon <= 0:
        raise ValueError("Positive tolerance required")
    terms = 1
    while True:
        if deadline is not None and perf_counter() > deadline:
            raise TimeoutError("Certification deadline reached")
        if max_terms is not None and terms > max_terms:
            raise TimeoutError("Certification series budget reached")
        interval = expression.interval(terms)
        if interval.width <= epsilon:
            return interval, terms
        terms *= 2


def clear_caches():
    natural_log_unit.cache_clear()
    log2_interval.cache_clear()
    natural_log_interval.cache_clear()
