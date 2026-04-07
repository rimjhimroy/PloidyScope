"""Compatibility helpers for legacy between-population analyses."""

from __future__ import annotations

from .diversity import calc_dxy_windows
from .rho import calc_rho_windows


def calc_bpm_windows(
    records, window_size=100000, minimum_snps=2, num_pops=2, outname="out"
):
    rho_rows, rho_summary = calc_rho_windows(
        records,
        window_size=window_size,
        minimum_snps=minimum_snps,
    )
    dxy_rows = {
        (row["pop1"], row["pop2"], row["chromosome"], row["window_pos_1"], row["window_pos_2"]): row
        for row in calc_dxy_windows(records, window_size=window_size)
    }

    legacy_rows = []
    for row in rho_rows:
        key = (
            row["pop1"],
            row["pop2"],
            row["chromosome"],
            row["window_pos_1"],
            row["window_pos_2"],
        )
        dxy_row = dxy_rows.get(key, {})
        legacy_rows.append(
            {
                "outname": outname,
                "scaff": row["chromosome"],
                "start": row["window_pos_1"],
                "end": row["window_pos_2"],
                "win_size": window_size,
                "num_sites": row["no_sites"],
                "num_snps": row["no_snps"],
                "Rho": row["rho"],
                "Fst": "NA",
                "dxy": dxy_row.get("avg_dxy", "NA"),
                "AFD": "NA",
                "FixedDiff": "NA",
            }
        )

    return legacy_rows, {"rho": rho_summary["pairs"]}
