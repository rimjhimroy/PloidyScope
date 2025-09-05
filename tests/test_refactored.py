import sys
from pathlib import Path

# Ensure test process can import the local package regardless of the directory name casing.
# The repository folder is named 'ScanTools_snakemake' but tests import 'scantools_snakemake'.
# Add the repository root to sys.path and, if needed, map the actual package module to
# the lowercase name so imports in tests work on case-sensitive filesystems.
# We need the directory that contains the package directory on sys.path. The tests
# live in <repo>/scripts/ScanTools_snakemake/tests; the package directory is
# <repo>/scripts/ScanTools_snakemake, so add its parent to sys.path.
root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root))
try:
    import scantools_snakemake  # preferred import
except ModuleNotFoundError:
    import importlib

    # On case-sensitive filesystems the actual directory is 'ScanTools_snakemake'.
    candidate = "ScanTools_snakemake"
    pkg_path = root / candidate
    if pkg_path.is_dir():
        mod = importlib.import_module(candidate)
        # Alias under the lowercase name tests expect.
        sys.modules["scantools_snakemake"] = mod

import tempfile
from scantools_snakemake.refactored.bpm import calc_bpm_windows
from scantools_snakemake.refactored.wpm import calc_wpm_windows


def make_synthetic_records():
    # create 2-population site with simple genotypes
    # format: pop, ploidy, scaff, pos, an, dp, genotypes...
    recs = []
    # site1 at pos 100: popA has 2 diploid individuals with acs 0 and 1; popB has 2 individuals acs 2,1
    recs.append(["POP1", "2", "chr1", "100", "4", "10", "0", "1"])
    recs.append(["POP2", "2", "chr1", "100", "4", "10", "2", "1"])
    # site2 at pos 200
    recs.append(["POP1", "2", "chr1", "200", "4", "10", "0", "0"])
    recs.append(["POP2", "2", "chr1", "200", "4", "10", "0", "1"])
    return recs


def test_bpm_basic():
    recs = make_synthetic_records()
    results, summary = calc_bpm_windows(
        recs, window_size=1000, minimum_snps=1, num_pops=2
    )
    assert isinstance(results, list)
    assert "rho" in summary


def test_wpm_basic():
    recs = make_synthetic_records()
    results = calc_wpm_windows(recs, sampind=2, window_size=1000, minimum_snps=1)
    assert isinstance(results, list)
    assert len(results) >= 1
    # Check keys and value types in first window
    win = results[0]
    expected_keys = [
        "pop",
        "ploidy",
        "sampind",
        "scaff",
        "start",
        "end",
        "window_size",
        "num_snps",
        "num_sites",
        "num_singletons",
        "avg_freq",
        "avg_Ehet",
        "Diversity",
        "ThetaW",
        "Pi",
        "ThetaH",
        "ThetaL",
        "D",
        "H",
        "E",
    ]
    for k in expected_keys:
        assert k in win
    assert isinstance(win["num_snps"], int)
    assert isinstance(win["Diversity"], float)
    # Edge case: all missing genotypes
    missing_recs = [
        ["POP1", "2", "chr1", "100", "4", "10", "-9", "-9"],
        ["POP1", "2", "chr1", "200", "4", "10", "-9", "-9"],
    ]
    results_missing = calc_wpm_windows(
        missing_recs, sampind=2, window_size=1000, minimum_snps=1
    )
    assert results_missing[0]["num_snps"] == 0
    assert results_missing[0]["Diversity"] == 0.0
