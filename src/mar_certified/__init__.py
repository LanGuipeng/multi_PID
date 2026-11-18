"""Certified Shannon MAR queries for finite discrete distributions."""

from .api import query
from .distributions import Distribution

__all__ = ["Distribution", "query"]
__version__ = "0.1.0"
