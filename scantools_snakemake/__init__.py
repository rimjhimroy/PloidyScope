"""scantools_snakemake package — canonical lowercase package entrypoint.

Expose the refactored API under `scantools_snakemake.refactored`. The
implementation lives in `refactored/` at the project root; this package
provides a stable import path used by tests and Snakemake wrappers.
"""

from . import refactored

__all__ = ["refactored"]

__version__ = "0.0.1"
