import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))

from ploidyscope.stats.diversity import calc_dxy_windows
from ploidyscope.stats.diversity import calc_pi_windows
from ploidyscope.stats.diversity import calc_tajima_d_windows
from ploidyscope.stats.rho import calc_rho_windows


def make_mixed_ploidy_records():
    return [
        ["POP1", "2", "chr1", "100", "4", "10", "0", "1"],
        ["POP2", "4", "chr1", "100", "8", "10", "1", "4"],
        ["POP1", "2", "chr1", "200", "4", "10", "1", "2"],
        ["POP2", "4", "chr1", "200", "8", "10", "0", "3"],
        ["POP1", "2", "chr1", "300", "4", "10", "0", "0"],
        ["POP2", "4", "chr1", "300", "8", "10", "0", "0"],
    ]


def test_rho_windows_mixed_ploidy():
    rows, summary = calc_rho_windows(
        make_mixed_ploidy_records(),
        window_size=1000,
        minimum_snps=1,
    )
    assert len(rows) == 1
    assert rows[0]["pop1"] == "POP1"
    assert rows[0]["pop2"] == "POP2"
    assert rows[0]["rho"] != "NA"
    assert rows[0]["rho_den"] != 0.0
    assert summary["pairs"][0]["rho"] != "NA"


def test_dxy_windows_arbitrary_ploidy_counts():
    rows = calc_dxy_windows(make_mixed_ploidy_records(), window_size=1000)
    assert len(rows) == 1
    assert rows[0]["count_comparisons"] > 0
    assert rows[0]["avg_dxy"] != "NA"


def test_pi_and_tajima_outputs_have_expected_fields():
    records = [
        ["POP1", "2", "chr1", "100", "4", "10", "0", "1", "2"],
        ["POP1", "2", "chr1", "200", "4", "10", "0", "0", "1"],
        ["POP1", "2", "chr1", "300", "4", "10", "0", "0", "0"],
    ]

    pi_rows = calc_pi_windows(records, window_size=1000)
    tajima_rows = calc_tajima_d_windows(records, window_size=1000)

    assert len(pi_rows) == 1
    assert len(tajima_rows) == 1
    assert "count_diffs" in pi_rows[0]
    assert "count_comparisons" in pi_rows[0]
    assert "raw_pi" in tajima_rows[0]
    assert "raw_watterson_theta" in tajima_rows[0]
    assert "tajima_d_stdev" in tajima_rows[0]
