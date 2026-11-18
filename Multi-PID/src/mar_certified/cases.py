from fractions import Fraction
from itertools import product
import random


def analytic_cases():
    result = {}

    def add(name, rows):
        table = {}
        for state, mass in rows:
            table[state] = table.get(state, Fraction()) + mass
        result[name] = table

    pairs = tuple(product(range(2), repeat=2))
    add("xor", [((a ^ b, a, b), Fraction(1, 4)) for a, b in pairs])
    add("copy", [((a, a, a), Fraction(1, 2)) for a in range(2)])
    add("unique", [((a, a, 0), Fraction(1, 2)) for a in range(2)])
    add("constant", [((0, a, b), Fraction(1, 4)) for a, b in pairs])
    add("independent", [((t, a, b), Fraction(1, 8)) for t, a, b in product(range(2), repeat=3)])
    add("and", [((a & b, a, b), Fraction(1, 4)) for a, b in pairs])
    add("or", [((a | b, a, b), Fraction(1, 4)) for a, b in pairs])
    add("double_copy", [((2 * a + b, a, b), Fraction(1, 4)) for a, b in pairs])
    add(
        "parity3",
        [((a ^ b ^ c, a, b, c), Fraction(1, 8)) for a, b, c in product(range(2), repeat=3)],
    )
    add("repeat3", [((a, a, a, a), Fraction(1, 2)) for a in range(2)])
    add("single", [((a, a), Fraction(1, 2)) for a in range(2)])
    rows = []
    for t, flip, erase in product(range(2), repeat=3):
        rows.append(
            (
                (t, t ^ flip, 2 if erase else t),
                Fraction(1, 2) * (Fraction(1, 10) if flip else Fraction(9, 10)) * Fraction(1, 2),
            )
        )
    add("bsc_bec", rows)
    for exponent in (4, 12, 30):
        e = Fraction(1, 2**exponent)
        add("rare_" + str(exponent), [((0, 0, 0), 1 - e), ((1, 1, 0), e)])
        rows = []
        for t, e1, e2 in product(range(2), repeat=3):
            a, b = Fraction(1, 10), Fraction(1, 10) + e / 10
            rows.append(
                ((t, t ^ e1, t ^ e2), Fraction(1, 2) * (a if e1 else 1 - a) * (b if e2 else 1 - b))
            )
        add("near_bsc_" + str(exponent), rows)
    add(
        "mixture",
        [((a ^ b, a, b), Fraction(1, 8)) for a, b in pairs]
        + [((2 + t, 2 + t, 2 + t), Fraction(1, 4)) for t in range(2)],
    )
    return result


def random_table(seed, n, target_size=3, source_size=2, sparsity=0.25):
    generator = random.Random(seed)
    counts = {
        state: (generator.randrange(1, 30) if generator.random() >= sparsity else 0)
        for state in product(range(target_size), *[range(source_size) for _ in range(n)])
    }
    if not sum(counts.values()):
        counts[next(iter(counts))] = 1
    total = sum(counts.values())
    return {state: Fraction(count, total) for state, count in counts.items()}


def input_record(table):
    return {
        "kind": "probabilities",
        "records": [{"state": list(s), "mass": str(p)} for s, p in sorted(table.items())],
    }
