import sys
from pathlib import Path

# Ensure local package import works when running this script directly
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))

from ploidyscope.stats.diversity import calc_tajima_d_windows
from ploidyscope.stats.rho import calc_rho_windows


def make_synthetic_records():
    recs = []
    recs.append(["POP1", "2", "chr1", "100", "4", "10", "0", "1"])
    recs.append(["POP2", "4", "chr1", "100", "8", "10", "1", "4"])
    recs.append(["POP1", "2", "chr1", "200", "4", "10", "1", "2"])
    recs.append(["POP2", "4", "chr1", "200", "8", "10", "0", "3"])
    return recs


if __name__ == "__main__":
    recs = make_synthetic_records()
    print("Calling calc_rho_windows...")
    results, summary = calc_rho_windows(recs, window_size=1000, minimum_snps=1)
    print("rho rows:", len(results), "summary pairs:", len(summary["pairs"]))
    print("Calling calc_tajima_d_windows...")
    wpm_res = calc_tajima_d_windows(recs[:2], window_size=1000)
    print("tajima rows:", len(wpm_res))
