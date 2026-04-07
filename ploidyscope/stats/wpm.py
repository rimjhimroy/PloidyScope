"""Compatibility helpers for legacy within-population analyses."""

from __future__ import annotations

from .diversity import calc_pi_windows
from .diversity import calc_tajima_d_windows


def calc_wpm_windows(records, sampind=5, window_size=50000, minimum_snps=2):
    pi_rows = {
        (row["pop"], row["chromosome"], row["window_pos_1"], row["window_pos_2"]): row
        for row in calc_pi_windows(records, window_size=window_size)
    }
    tajima_rows = {
        (row["pop"], row["chromosome"], row["window_pos_1"], row["window_pos_2"]): row
        for row in calc_tajima_d_windows(records, window_size=window_size)
    }

    results = []
    for key in sorted(pi_rows):
        pi_row = pi_rows[key]
        tajima_row = tajima_rows.get(key, {})
        results.append(
            {
                "pop": pi_row["pop"],
                "scaff": pi_row["chromosome"],
                "start": pi_row["window_pos_1"],
                "end": pi_row["window_pos_2"],
                "window_size": window_size,
                "num_sites": pi_row["no_sites"],
                "Diversity": pi_row["pi"],
                "Pi": tajima_row.get("raw_pi", "NA"),
                "ThetaW": tajima_row.get("raw_watterson_theta", "NA"),
                "D": tajima_row.get("tajima_d", "NA"),
            }
        )
    return results
