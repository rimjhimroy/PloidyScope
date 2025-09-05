import sys
from pathlib import Path

# Ensure local package import works when running this script directly
root = Path(__file__).resolve().parents[2]
# Add the directory that contains the package folder (i.e. the 'scripts' directory)
sys.path.insert(0, str(root))
try:
    import scantools_snakemake
except ModuleNotFoundError:
    import importlib

    candidate = "ScanTools_snakemake"
    pkg_path = root / candidate
    if pkg_path.is_dir():
        mod = importlib.import_module(candidate)
        sys.modules["scantools_snakemake"] = mod

from scantools_snakemake.refactored.bpm import calc_bpm_windows
from scantools_snakemake.refactored.wpm import calc_wpm_windows


def make_synthetic_records():
    recs = []
    recs.append(["POP1", "2", "chr1", "100", "4", "10", "0", "1"])
    recs.append(["POP2", "2", "chr1", "100", "4", "10", "2", "1"])
    recs.append(["POP1", "2", "chr1", "200", "4", "10", "0", "0"])
    recs.append(["POP2", "2", "chr1", "200", "4", "10", "0", "1"])
    return recs


if __name__ == "__main__":
    recs = make_synthetic_records()
    print("Calling calc_bpm_windows...")
    results, summary = calc_bpm_windows(
        recs, window_size=1000, minimum_snps=1, num_pops=2
    )
    print(
        "bpm results type:",
        type(results),
        "summary keys:",
        list(summary.keys()) if isinstance(summary, dict) else summary,
    )
    print("Calling calc_wpm_windows...")
    wpm_res = calc_wpm_windows(recs, sampind=2, window_size=1000, minimum_snps=1)
    print("wpm results type:", type(wpm_res), "len:", len(wpm_res))
