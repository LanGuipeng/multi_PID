from itertools import combinations
import json


def contained(first, second):
    return first & second == first


def validate_antichain(alpha, n_sources):
    alpha = tuple(alpha)
    if type(n_sources) is not int or n_sources < 1:
        raise ValueError("Positive source count required")
    if not alpha or len(set(alpha)) != len(alpha):
        raise ValueError("Nonempty antichain without duplicate sets required")
    if any(type(mask) is not int or not 0 < mask < 1 << n_sources for mask in alpha):
        raise ValueError("Source mask outside universe")
    if any(contained(a, b) or contained(b, a) for a, b in combinations(alpha, 2)):
        raise ValueError("Comparable sets do not form an antichain")
    return tuple(sorted(alpha))


def parse_antichain(value, n_sources):
    if isinstance(value, str):
        if value.lstrip().startswith("["):
            value = json.loads(value)
        else:
            value = [
                [int(token.strip()) for token in group.split(",")] for group in value.split("|")
            ]
    if isinstance(value, tuple) and all(type(mask) is int for mask in value):
        return validate_antichain(value, n_sources)
    if not isinstance(value, (list, tuple)):
        raise ValueError("Antichain must be string, source-list JSON or mask tuple")
    masks = []
    for group in value:
        if not isinstance(group, (list, tuple)) or not group or len(set(group)) != len(group):
            raise ValueError("Each source set is a nonempty list of distinct 1-based indices")
        if any(type(i) is not int or not 1 <= i <= n_sources for i in group):
            raise ValueError("Source index outside 1..n")
        masks.append(sum(1 << (i - 1) for i in group))
    return validate_antichain(masks, n_sources)


def boundary_sets(alpha, n_sources, max_boundary_subsets=1048576, check=None):
    alpha = validate_antichain(alpha, n_sources)
    universe = (1 << n_sources) - 1
    if alpha == tuple(1 << i for i in range(n_sources)):
        return ()
    if len(alpha) == 1 and alpha[0].bit_count() == 1:
        return (universe ^ alpha[0],) if n_sources > 1 else ()
    if alpha == (universe,):
        return tuple(sorted(universe ^ (1 << i) for i in range(n_sources))) if n_sources > 1 else ()
    if universe > max_boundary_subsets:
        raise TimeoutError("Boundary subset budget exceeded")
    outside = set()
    for mask in range(1, universe + 1):
        if check and mask % 256 == 0:
            check()
        if not any(contained(a, mask) for a in alpha):
            outside.add(mask)
    return tuple(
        sorted(
            mask
            for mask in outside
            if all(mask | (1 << i) not in outside for i in range(n_sources) if not mask & (1 << i))
        )
    )


def leq(alpha, beta):
    return all(any(contained(a, b) for a in alpha) for b in beta)


def antichains(n_sources, maximum_sources=3):
    if not 1 <= n_sources <= maximum_sources or maximum_sources > 4:
        raise ValueError("Full lattice allowed only for n <= configured maximum <= 4")

    def extend(prefix, candidates):
        for index, mask in enumerate(candidates):
            current = prefix + (mask,)
            yield current
            remaining = tuple(
                other
                for other in candidates[index + 1 :]
                if not contained(mask, other) and not contained(other, mask)
            )
            yield from extend(current, remaining)

    return tuple(extend((), tuple(range(1, 1 << n_sources))))


def topological_antichains(n_sources, maximum_sources=3):
    pending = set(antichains(n_sources, maximum_sources))
    result = []
    while pending:
        ready = sorted(
            beta
            for beta in pending
            if all(gamma == beta or not leq(gamma, beta) for gamma in pending)
        )
        if not ready:
            raise ArithmeticError("Order cycle")
        result.extend(ready)
        pending.difference_update(ready)
    return tuple(result)
