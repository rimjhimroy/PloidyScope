#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

# ensure local ScanTools package is importable
pkg_dir = Path(__file__).resolve().parents[0] / "ScanTools"
sys.path.insert(0, str(pkg_dir))

# The original ScanTools package contains modules like ScanTools.py, bpm.py, wpm.py, etc.
# This wrapper selects a function to run based on `--mode` and passes through arguments.


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["bpm", "wpm", "calcAFS", "calcFreqs", "MissingData", "repol"],
        default="bpm",
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--popkey", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--window", required=True, help="chr:start-end")
    args = parser.parse_args()

    mode = args.mode
    vcf = args.vcf
    popkey = args.popkey
    out = args.out
    window = args.window

    # Import the relevant module dynamically
    if mode == "bpm":
        import bpm as mod
    elif mode == "wpm":
        import wpm as mod
    elif mode == "calcAFS":
        import calcAFS as mod
    elif mode == "calcFreqs":
        import calcFreqs_atSites as mod
    elif mode == "MissingData":
        import MissingData as mod
    elif mode == "repol":
        import repol as mod
    else:
        raise SystemExit("Unknown mode")

    # Many ScanTools scripts expect command-line usage; try to call main() if present
    if hasattr(mod, "main"):
        argv_backup = sys.argv
        sys.argv = [
            argv_backup[0],
            "--vcf",
            vcf,
            "--popkey",
            popkey,
            "--out",
            out,
            "--window",
            window,
        ]
        try:
            mod.main()
        finally:
            sys.argv = argv_backup
    else:
        # fallback: call a function if present
        if hasattr(mod, "run"):
            mod.run(vcf=vcf, popkey=popkey, out=out, window=window)
        else:
            raise SystemExit("Module has no callable entrypoint")


if __name__ == "__main__":
    main()
