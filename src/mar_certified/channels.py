from fractions import Fraction


def validate_channel(columns):
    columns = tuple((Fraction(a), Fraction(b)) for a, b in columns)
    if any(a < 0 or b < 0 for a, b in columns):
        raise ValueError("Negative channel probability")
    if sum(a for a, b in columns) != 1 or sum(b for a, b in columns) != 1:
        raise ValueError("Each conditional channel row must sum to one")
    return tuple((a, b) for a, b in columns if a or b)


def decision_value(columns, lam):
    lam = Fraction(lam)
    if lam < 0:
        raise ValueError("Nonnegative decision parameter required")
    return sum((max(Fraction(), a - lam * b) for a, b in columns), Fraction()) - max(
        Fraction(), 1 - lam
    )


def affine_line(columns, probe):
    active = tuple((a, b) for a, b in columns if a - probe * b > 0)
    return (
        sum((a for a, b in active), Fraction()) - (1 if probe < 1 else 0),
        -sum((b for a, b in active), Fraction()) + (1 if probe < 1 else 0),
    )


def mutual_information_expression(distribution, mask):
    from .intervals import LogExpression

    coordinates = distribution.coordinates(mask)
    joint = distribution.marginal((0,) + coordinates)
    source = distribution.marginal(coordinates)
    priors = distribution.target_probabilities()
    return LogExpression.from_terms(
        (mass, mass / (priors[state[0]] * source[state[1:]]))
        for state, mass in joint.items()
        if mass
    )
