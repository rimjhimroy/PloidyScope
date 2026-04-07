#!/usr/bin/env python3
"""Run PloidyScope statistics on a recoded table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .diversity import calc_dxy_windows
from .diversity import calc_pi_windows
from .diversity import calc_tajima_d_windows
from .rho import calc_rho_windows


def load_records(path: str) -> list[list[str]]:
    records = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line:
                records.append(line.split("\t"))
    return records


def write_tsv(rows: list[dict[str, object]], path: str) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        output_path.write_text("", encoding="utf-8")
        return

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stat",
        choices=["rho", "dxy", "pi", "tajima_d"],
        required=True,
    )
    parser.add_argument("--infile", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--window-size", type=int, default=100000)
    parser.add_argument("--minimum-snps", type=int, default=2)
    parser.add_argument("--populations", nargs="*", default=None)
    args = parser.parse_args()

    records = load_records(args.infile)

    if args.stat == "rho":
        rows, _ = calc_rho_windows(
            records,
            window_size=args.window_size,
            minimum_snps=args.minimum_snps,
            populations=args.populations,
        )
    elif args.stat == "dxy":
        rows = calc_dxy_windows(
            records,
            window_size=args.window_size,
            populations=args.populations,
        )
    elif args.stat == "pi":
        rows = calc_pi_windows(
            records,
            window_size=args.window_size,
            populations=args.populations,
        )
    else:
        rows = calc_tajima_d_windows(
            records,
            window_size=args.window_size,
            populations=args.populations,
        )

    write_tsv(rows, args.out)


if __name__ == "__main__":
    main()
