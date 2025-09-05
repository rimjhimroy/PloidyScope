"""Refactored ScanTools API (canonical package path)

This module exposes the functions implemented in the project's top-level
`refactored` package so `scantools_snakemake.refactored` works as an import
path without duplicating code.
"""

# Import from the top-level `refactored` package that lives at
# scripts/ScanTools_snakemake/refactored
from refactored.bpm import calc_bpm_windows
from refactored.wpm import calc_wpm_windows

__all__ = ["calc_bpm_windows", "calc_wpm_windows"]
