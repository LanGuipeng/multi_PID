from fractions import Fraction
from itertools import combinations
from .channels import validate_channel, affine_line


def upper_envelope(lines, left, right):
    unique = []
    for intercept, slope in sorted(lines, key=lambda line: (line[1], line[0])):
        if unique and unique[-1][1] == slope:
            unique[-1] = (intercept, slope)
        else:
            unique.append((intercept, slope))
    hull = []
    for intercept, slope in unique:
        start = None
        while hull:
            previous_intercept, previous_slope, previous_start = hull[-1]
            start = (previous_intercept - intercept) / (slope - previous_slope)
            if previous_start is None or start > previous_start:
                break
            hull.pop()
        if not hull:
            start = None
        hull.append((intercept, slope, start))
    segments = []
    for index, (intercept, slope, start) in enumerate(hull):
        stop = hull[index + 1][2] if index + 1 < len(hull) else None
        clipped_left = max(left, start) if start is not None else left
        clipped_right = stop if right is None else right if stop is None else min(right, stop)
        if clipped_right is None or clipped_left < clipped_right:
            segments.append((clipped_left, clipped_right, intercept, slope))
    return tuple(segments)


def lower_envelope(lines, left, right):
    return tuple(
        (a, b, -c, -d)
        for a, b, c, d in upper_envelope(tuple((-c, -d) for c, d in lines), left, right)
    )


def append_segment(result, segment):
    left, right, intercept, slope = segment
    if intercept == 0 and slope == 0:
        return
    if result and result[-1][1] == left and result[-1][2:] == (intercept, slope):
        result[-1] = (result[-1][0], right, intercept, slope)
    else:
        result.append(segment)


def positive_difference(primary, competitor):
    result = []
    first = second = 0
    while first < len(primary) and second < len(competitor):
        a, b, c, d = primary[first]
        e, f, g, h = competitor[second]
        left = max(a, e)
        right = f if b is None else b if f is None else min(b, f)
        intercept, slope = c - g, d - h
        points = [left]
        if slope:
            crossing = -intercept / slope
            if crossing > left and (right is None or crossing < right):
                points.append(crossing)
        if right is not None:
            points.append(right)
        for index, start in enumerate(points):
            end = points[index + 1] if index + 1 < len(points) else right
            if end == start or (end is None and right is not None):
                continue
            sample = (start + end) / 2 if end is not None else start + 1
            if intercept + slope * sample > 0:
                append_segment(result, (start, end, intercept, slope))
        if right is None:
            break
        if b == right:
            first += 1
        if f == right:
            second += 1
    return tuple(result)


def optimized_segments(channels, primary_count, atom=True, check=None, max_segments=100000):
    if not channels or not 1 <= primary_count <= len(channels):
        raise ValueError("Nonempty primary family required")
    events = [(Fraction(1), None, Fraction(), Fraction())]
    lines = []
    for index, channel in enumerate(channels):
        if sum(a for a, b in channel) != 1 or sum(b for a, b in channel) != 1:
            raise ValueError("Channel rows must sum to one")
        intercept = sum((a for a, b in channel if a > 0), Fraction()) - 1
        slope = 1 - sum((b for a, b in channel if a > 0), Fraction())
        lines.append((intercept, slope))
        for a, b in channel:
            if a < 0 or b < 0:
                raise ValueError("Negative channel entry")
            if a and b:
                events.append((a / b, index, a, b))
    groups = [(Fraction(), [])]
    for point, index, a, b in sorted(events, key=lambda event: event[0]):
        if groups[-1][0] != point:
            groups.append((point, []))
        if index is not None:
            groups[-1][1].append((index, a, b))
    knots = [point for point, group in groups]
    result = []
    for position, left in enumerate(knots):
        if check:
            check()
        if len(result) > max_segments:
            raise TimeoutError("Segment budget exceeded")
        if left:
            for index, a, b in groups[position][1]:
                intercept, slope = lines[index]
                lines[index] = intercept - a, slope + b
            if left == 1:
                lines = [(intercept + 1, slope - 1) for intercept, slope in lines]
        right = knots[position + 1] if position + 1 < len(knots) else None
        primary = lower_envelope(lines[:primary_count], left, right)
        if not atom:
            for segment in primary:
                append_segment(result, segment)
            continue
        competitor = (
            upper_envelope(lines[primary_count:], left, right)
            if primary_count < len(lines)
            else ((left, right, Fraction(), Fraction()),)
        )
        for segment in positive_difference(primary, competitor):
            append_segment(result, segment)
    return tuple(result)


def enumerated_segments(channels, primary_count, atom=True, check=None, max_segments=100000):
    channels = tuple(validate_channel(channel) for channel in channels)
    if not channels or not 1 <= primary_count <= len(channels):
        raise ValueError("Invalid primary channel count")
    knots = sorted(
        {Fraction(), Fraction(1)} | {a / b for channel in channels for a, b in channel if a and b}
    )
    result = []
    for position, left in enumerate(knots):
        if check:
            check()
        right = knots[position + 1] if position + 1 < len(knots) else None
        probe = (left + right) / 2 if right is not None else left + 1
        lines = tuple(affine_line(channel, probe) for channel in channels)
        if right is None and any(d for c, d in lines):
            raise ArithmeticError("Finite channel final tail must be constant")
        points = {left}
        if right is not None:
            points.add(right)
        for count, (first, second) in enumerate(
            combinations(lines + ((Fraction(), Fraction()),), 2)
        ):
            if check and count % 256 == 0:
                check()
            if first[1] != second[1]:
                crossing = (second[0] - first[0]) / (first[1] - second[1])
                if crossing > left and (right is None or crossing < right):
                    points.add(crossing)
        points = sorted(points)
        for index, start in enumerate(points):
            end = points[index + 1] if index + 1 < len(points) else right
            if end == start:
                continue
            sample = (start + end) / 2 if end is not None else start + 1
            primary = min(lines[:primary_count], key=lambda line: line[0] + sample * line[1])
            competitor = (
                max(lines[primary_count:], key=lambda line: line[0] + sample * line[1])
                if atom and primary_count < len(lines)
                else (Fraction(), Fraction())
            )
            c, d = primary[0] - competitor[0], primary[1] - competitor[1]
            if c + d * sample > 0:
                append_segment(result, (start, end, c, d))
            if len(result) > max_segments:
                raise TimeoutError("Segment budget exceeded")
    return tuple(result)


def decision_segments(
    channels, primary_count, atom=True, algorithm="enumerate", check=None, max_segments=100000
):
    if algorithm == "enumerate":
        return enumerated_segments(channels, primary_count, atom, check, max_segments)
    if algorithm == "optimized":
        return optimized_segments(channels, primary_count, atom, check, max_segments)
    raise ValueError("Algorithm must be enumerate or optimized")
