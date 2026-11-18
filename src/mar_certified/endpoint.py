from fractions import Fraction
from time import perf_counter
from .intervals import Interval, log2_interval, weight_interval
from .lattice import boundary_sets, validate_antichain


def curve_value(channel, value):
    value = Fraction(value)
    return sum((max(Fraction(), a - value * b) for a, b in channel), Fraction()) - max(
        Fraction(), 1 - value
    )


def endpoint_bounds(channels, source_count, left, right, kind="atom"):
    left, right = Fraction(left), Fraction(right)
    if not 0 <= left < right or left < 1 < right:
        raise ValueError("A cell must have ordered nonnegative endpoints and cannot cross one")
    bounds = tuple(
        (
            min(curve_value(channel, left), curve_value(channel, right)),
            max(curve_value(channel, left), curve_value(channel, right)),
        )
        for channel in channels
    )
    if not 1 <= source_count <= len(bounds):
        raise ValueError("Invalid source count")
    lower = min(item[0] for item in bounds[:source_count])
    upper = min(item[1] for item in bounds[:source_count])
    if kind == "atom":
        lower = max(
            Fraction(), lower - max((item[1] for item in bounds[source_count:]), default=Fraction())
        )
        upper = max(
            Fraction(), upper - max((item[0] for item in bounds[source_count:]), default=Fraction())
        )
    elif kind != "cumulative":
        raise ValueError("Unknown query kind")
    return lower, upper


def dyadic_cells(level):
    if not isinstance(level, int) or level < 1:
        raise ValueError("Positive grid level required")
    subdivisions = 1 << level
    for exponent in range(-level, level):
        start = Fraction(2) ** exponent
        step = start / subdivisions
        for index in range(subdivisions):
            yield start + index * step, start + (index + 1) * step


def tail_bound(priors, left, right):
    priors = tuple(Fraction(p) for p in priors)
    if any(not 0 < p < 1 for p in priors):
        raise ValueError("Tail bound requires nonconstant positive priors")
    left, right = Fraction(left), Fraction(right)
    if not 0 < left <= 1 <= right:
        raise ValueError("Tail endpoints must straddle one")
    return Fraction(3, 2) * left + Fraction(3, 4) / right**2 * sum(
        ((1 - p) ** 2 / p for p in priors), Fraction()
    )


def floor_scaled(value, scale):
    value = Fraction(value) * scale
    return value.numerator // value.denominator


def ceil_scaled(value, scale):
    value = Fraction(value) * scale
    return -((-value.numerator) // value.denominator)


def evaluate(
    distribution,
    alpha,
    kind="atom",
    tol="1e-8",
    max_seconds=10.0,
    max_cells=4096,
    max_level=20,
    max_terms=256,
    **options,
):
    started = perf_counter()
    epsilon = Fraction(tol)
    if epsilon <= 0:
        raise ValueError("Positive tolerance required")
    if kind not in ("atom", "cumulative"):
        raise ValueError("Unknown query kind")
    if max_seconds is not None and max_seconds < 0:
        raise ValueError("Time budget must be nonnegative")
    if not isinstance(max_cells, int) or max_cells < 0 or max_level < 0 or max_terms < 1:
        raise ValueError("Invalid grid or series budget")
    alpha = validate_antichain(alpha, distribution.n_sources)
    priors = {
        target: prior for target, prior in distribution.target_probabilities().items() if prior > 0
    }
    timings = dict(
        boundary=0.0, marginalization=0.0, envelope=0.0, integration=0.0, certification=0.0
    )
    metrics = dict(
        n_sources=distribution.n_sources,
        support_size=len(distribution.table),
        alpha_size=len(alpha),
        boundary_size=0,
        channel_columns=0,
        grid_cells=0,
        completed_level=0,
        series_terms=0,
        evaluated_cells=0,
        method="global_dyadic_endpoint",
    )
    if len(priors) == 1:
        return dict(
            expression=None,
            interval=Interval.point(0),
            tolerance_met=True,
            exact_zero=True,
            strictly_positive=False,
            timings={**timings, "total": perf_counter() - started},
            metrics=metrics,
            status="success",
            guarantee="deterministic_rational_enclosure",
            error_budget=dict(
                tail=Fraction(),
                endpoint=Fraction(),
                weight_arithmetic=Fraction(),
                summation_rounding=Fraction(),
                total_width=Fraction(),
            ),
        )
    deadline = float("inf") if max_seconds is None else started + max_seconds
    upper_global = log2_interval(Fraction(len(priors)), 8).upper
    best = Interval(0, upper_global)

    def check_deadline():
        if perf_counter() >= deadline:
            raise TimeoutError("time_budget")

    before = perf_counter()
    try:
        check_deadline()
        competitors = (
            boundary_sets(
                alpha,
                distribution.n_sources,
                options.get("max_boundary_subsets", 1048576),
                check_deadline,
            )
            if kind == "atom"
            else ()
        )
        timings["boundary"] = perf_counter() - before
        metrics["boundary_size"] = len(competitors)
        masks = alpha + competitors
        before = perf_counter()
        channels_by_target = {}
        for target in priors:
            check_deadline()
            channels_by_target[target] = tuple(distribution.channel(mask, target) for mask in masks)
        timings["marginalization"] = perf_counter() - before
        metrics["channel_columns"] = sum(
            len(channel) for channels in channels_by_target.values() for channel in channels
        )
    except TimeoutError as error:
        return dict(
            expression=None,
            interval=best,
            tolerance_met=best.width <= epsilon,
            exact_zero=False,
            strictly_positive=False,
            timings={**timings, "total": perf_counter() - started},
            metrics=metrics,
            status="success" if best.width <= epsilon else "timeout",
            reason=str(error),
            guarantee="deterministic_rational_enclosure",
            error_budget=dict(
                tail=None,
                endpoint=None,
                weight_arithmetic=None,
                summation_rounding=None,
                total_width=best.width,
            ),
            history=[],
        )
    errors = dict(
        tail=None,
        endpoint=None,
        weight_arithmetic=None,
        summation_rounding=None,
        total_width=best.width,
    )
    history = []
    reason = "maximum_level"
    for level in range(1, max_level + 1):
        count = 2 * level * (1 << level)
        if count > max_cells:
            reason = "grid_budget"
            break
        if perf_counter() >= deadline:
            reason = "time_budget"
            break
        left_tail, right_tail = Fraction(1, 1 << level), Fraction(1 << level)
        tail = tail_bound(priors.values(), left_tail, right_tail)
        weight_tolerance = Fraction(1, (1 << level) * len(priors) * count)
        rounding_bits = level + (len(priors) * count).bit_length() + 20
        scale = 1 << rounding_bits
        lower_sum = upper_sum = enclosure_sum = weight_error_sum = 0
        maximum_terms = 0
        interrupted = False
        for left, right in dyadic_cells(level):
            for target, prior in priors.items():
                if perf_counter() >= deadline:
                    reason = "time_budget"
                    interrupted = True
                    break
                before = perf_counter()
                lower, upper = endpoint_bounds(
                    channels_by_target[target], len(alpha), left, right, kind
                )
                timings["envelope"] += perf_counter() - before
                metrics["evaluated_cells"] += 1
                if upper == 0:
                    continue
                terms = 1
                before = perf_counter()
                while True:
                    if terms > max_terms:
                        reason = "series_budget"
                        interrupted = True
                        break
                    if perf_counter() >= deadline:
                        reason = "time_budget"
                        interrupted = True
                        break
                    weight = weight_interval(prior, left, right, terms)
                    if weight.width <= weight_tolerance:
                        break
                    terms *= 2
                timings["certification"] += perf_counter() - before
                if interrupted:
                    break
                maximum_terms = max(maximum_terms, terms)
                before = perf_counter()
                lower_sum += floor_scaled(prior * lower * weight.lower, scale)
                upper_sum += ceil_scaled(prior * upper * weight.upper, scale)
                enclosure_sum += ceil_scaled(prior * (upper - lower) * weight.lower, scale)
                weight_error_sum += ceil_scaled(prior * upper * weight.width, scale)
                timings["integration"] += perf_counter() - before
            if interrupted:
                break
        if interrupted:
            break
        round_interval = Interval(Fraction(lower_sum, scale), Fraction(upper_sum, scale) + tail)
        best = Interval(
            max(best.lower, round_interval.lower), min(best.upper, round_interval.upper)
        )
        errors = dict(
            tail=tail,
            endpoint=Fraction(enclosure_sum, scale),
            weight_arithmetic=Fraction(weight_error_sum, scale),
            summation_rounding=Fraction(2 * count * len(priors), scale),
            total_width=best.width,
        )
        metrics.update(
            grid_cells=count,
            completed_level=level,
            series_terms=maximum_terms,
            rounding_bits=rounding_bits,
        )
        history.append(
            dict(level=level, grid_cells=count, interval=best, error_budget=errors.copy())
        )
        if best.width <= epsilon:
            reason = None
            break
    timings["total"] = perf_counter() - started
    metrics["cache_hits"] = distribution.cache_hits
    metrics["cache_misses"] = distribution.cache_misses
    success = best.width <= epsilon
    return dict(
        expression=None,
        interval=best,
        tolerance_met=success,
        exact_zero=False,
        strictly_positive=best.lower > 0,
        timings=timings,
        metrics=metrics,
        status="success" if success else "timeout",
        reason=reason,
        guarantee="deterministic_rational_enclosure",
        error_budget=errors,
        history=history,
        series_terms=metrics["series_terms"],
    )
