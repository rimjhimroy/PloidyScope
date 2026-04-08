#!/usr/bin/env python3
"""Run PloidyScope statistics on a recoded table or VCF."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .common import iter_loci_from_vcf
from .diversity import calc_dxy_windows
from .diversity import calc_hudson_fst_windows
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


def load_input_records(
    infile: str | None,
    vcf: str | None,
    popmap: str | None,
    populations: list[str] | None,
):
    if vcf:
        if not popmap:
            raise ValueError("--popmap is required when using --vcf")
        return iter_loci_from_vcf(vcf, popmap, populations=populations)
    if infile:
        return load_records(infile)
    raise ValueError("Provide either --infile or --vcf")


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
        choices=["rho", "dxy", "pi", "tajima_d", "hudson"],
        required=True,
    )
    parser.add_argument("--infile")
    parser.add_argument("--vcf")
    parser.add_argument("--popmap")
    parser.add_argument("--out", required=True)
    parser.add_argument("--window-size", type=int, default=100000)
    parser.add_argument("--minimum-snps", type=int, default=2)
    parser.add_argument("--populations", nargs="*", default=None)
    args = parser.parse_args()

    if bool(args.infile) == bool(args.vcf):
        raise SystemExit("Provide exactly one of --infile or --vcf")

    records = load_input_records(args.infile, args.vcf, args.popmap, args.populations)

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
    elif args.stat == "hudson":
        rows = calc_hudson_fst_windows(
            records,
            window_size=args.window_size,
            populations=args.populations,
            minimum_snps=args.minimum_snps,
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
