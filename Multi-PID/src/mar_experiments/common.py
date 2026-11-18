"""Shared output handling for serial reproducibility runs."""

import json
import sys
from functools import wraps
from pathlib import Path


def decimal_output(function):
    """Allow large exact certificates during a run, then restore Python's limit."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        previous = None
        if hasattr(sys, "get_int_max_str_digits"):
            previous = sys.get_int_max_str_digits()
            sys.set_int_max_str_digits(0)
        try:
            return function(*args, **kwargs)
        finally:
            if previous is not None:
                sys.set_int_max_str_digits(previous)

    return wrapped


def write_json(output_dir, name, data):
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(data, indent=2, default=str, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)
    return path
