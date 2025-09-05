"""Refactored ScanTools API
Expose pure functions for core analyses (BPM, WPM) that accept iterables or pandas DataFrames and return Python dicts/lists.
"""

from .bpm import calc_bpm_windows
from .wpm import calc_wpm_windows

__all__ = ["calc_bpm_windows", "calc_wpm_windows"]
