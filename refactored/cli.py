#!/usr/bin/env python3
"""CLI to run refactored ScanTools functions on a tabular input file.
This is a minimal runner: it accepts an input table with same columns the original scripts expected.
"""
import argparse
from pathlib import Path
from .bpm import calc_bpm_windows
from .wpm import calc_wpm_windows

parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["bpm", "wpm"], default="bpm")
parser.add_argument("--infile", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--window", type=int, default=100000)
args = parser.parse_args()

records = []
with open(args.infile) as fh:
    for line in fh:
        line = line.strip().split("\t")
        if not line:
            continue
        records.append(line)

if args.mode == "bpm":
    res = calc_bpm_windows(records, window_size=args.window)
else:
    res = calc_wpm_windows(records, window_size=args.window)

# write simple TSV
outp = Path(args.out)
with outp.open("w") as oh:
    # naive header
    oh.write("scaff\tstart\tend\twin_size\tnum_sites\n")
    for i, r in enumerate(res):
        oh.write(f'chr\t0\t0\t{args.window}\t{r.get("num_sites",0)}\n')
print("Wrote", outp)
