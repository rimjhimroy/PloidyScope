import sys
from pathlib import Path
import textwrap

import pytest

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))

from ploidyscope.stats.common import build_record
from ploidyscope.stats.common import iter_loci_from_vcf
from ploidyscope.stats.diversity import calc_dxy_windows
from ploidyscope.stats.diversity import calc_dxy_windows_scantools_mean
from ploidyscope.stats.diversity import calc_hudson_fst_windows
from ploidyscope.stats.diversity import calc_pi_windows
from ploidyscope.stats.diversity import calc_tajima_d_windows
from ploidyscope.stats.rho import calc_diploid_wc_fst_site
from ploidyscope.stats.rho import calc_rho_windows
from scripts.compare_real_fixture import build_comparison_plan
from scripts.compare_real_fixture import write_svg_heatmap
from scripts.compare_real_fixture import normalize_scantools_bpm_rows
from scripts.compare_real_fixture import parse_metadata
from scripts.compare_real_fixture import record_passes_ploidyscope_vcf_filters
from scripts.compare_real_fixture import SampleMetadata


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


def test_diploid_wc_fst_site_keeps_c_only_sparse_site_contribution():
    sparse_site = [
        build_record("POP1", 2, "chr1", 100, [None, None]),
        build_record("POP2", 2, "chr1", 100, [2, 1]),
    ]
    shared_site = [
        build_record("POP1", 2, "chr1", 200, [0, 0]),
        build_record("POP2", 2, "chr1", 200, [1, 1]),
    ]

    sparse_components = calc_diploid_wc_fst_site(sparse_site)
    shared_components = calc_diploid_wc_fst_site(shared_site)
    rows, summary = calc_rho_windows(sparse_site + shared_site, window_size=1000, minimum_snps=1)

    assert sparse_components is not None
    assert shared_components is not None
    assert sparse_components.polymorphic is True
    assert sparse_components.fst_num == pytest.approx(0.0)
    assert sparse_components.fst_den == pytest.approx(0.5)
    assert shared_components.fst_num == pytest.approx(0.25)
    assert shared_components.fst_den == pytest.approx(0.5)
    assert rows[0]["fst_no_snps"] == 2
    assert summary["pairs"][0]["fst_no_snps"] == 2
    assert summary["pairs"][0]["fst"] == pytest.approx(0.25)


def test_dxy_windows_arbitrary_ploidy_counts():
    rows = calc_dxy_windows(make_mixed_ploidy_records(), window_size=1000)
    assert len(rows) == 1
    assert rows[0]["count_comparisons"] > 0
    assert rows[0]["avg_dxy"] != "NA"


def test_dxy_scantools_mean_uses_unweighted_site_average():
    records = [
        ["POP1", "2", "chr1", "100", "4", "10", "0", "2"],
        ["POP2", "4", "chr1", "100", "8", "10", "0", "4"],
        ["POP1", "2", "chr1", "200", "4", "10", "0", "1"],
        ["POP2", "4", "chr1", "200", "8", "10", "4", "-9"],
    ]

    weighted = calc_dxy_windows(records, window_size=1000)
    unweighted = calc_dxy_windows_scantools_mean(records, window_size=1000)

    assert len(weighted) == 1
    assert len(unweighted) == 1
    assert weighted[0]["avg_dxy"] == pytest.approx(0.5833333333333334)
    assert unweighted[0]["dxy_scantools_mean"] == pytest.approx(0.625)
    assert weighted[0]["avg_dxy"] != unweighted[0]["dxy_scantools_mean"]


def test_hudson_fst_windows_has_expected_value():
    rows = calc_hudson_fst_windows(make_mixed_ploidy_records(), window_size=1000, minimum_snps=1)

    assert len(rows) == 1
    assert rows[0]["avg_hudson_fst"] == pytest.approx(0.07936507936507936)
    assert rows[0]["no_snps"] == 2


def test_hudson_fst_windows_supports_tetraploid_pairs():
    records = [
        ["POP1", "4", "chr1", "100", "8", "10", "0", "2"],
        ["POP2", "4", "chr1", "100", "8", "10", "1", "4"],
        ["POP1", "4", "chr1", "200", "8", "10", "0", "0"],
        ["POP2", "4", "chr1", "200", "8", "10", "4", "4"],
    ]

    rows = calc_hudson_fst_windows(records, window_size=1000, minimum_snps=1)

    assert len(rows) == 1
    assert rows[0]["avg_hudson_fst"] != "NA"
    assert rows[0]["no_snps"] == 2


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


def test_iter_loci_from_vcf_groups_samples_by_population(tmp_path):
    vcf_path = tmp_path / "toy.vcf"
    popmap_path = tmp_path / "populations.tsv"

    vcf_path.write_text(
        textwrap.dedent(
            """\
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
            chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1/1/1\t0/0/0/1
            chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0/0/0\t0/1/1/1
            """
        ),
        encoding="utf-8",
    )
    popmap_path.write_text("sample\tpopulation\ns1\tPOP1\ns2\tPOP1\ns3\tPOP2\ns4\tPOP2\n", encoding="utf-8")

    loci = list(iter_loci_from_vcf(str(vcf_path), str(popmap_path)))

    assert len(loci) == 2
    assert [record.population for record in loci[0]] == ["POP1", "POP2"]
    assert loci[0][0].ploidy == 2
    assert loci[0][1].ploidy == 4
    assert loci[0][0].observed_alt == 1
    assert loci[0][1].observed_alt == 4


def test_vcf_records_run_through_rho_windows(tmp_path):
    vcf_path = tmp_path / "toy.vcf"
    popmap_path = tmp_path / "populations.tsv"

    vcf_path.write_text(
        textwrap.dedent(
            """\
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
            chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1/1/1\t0/0/0/1
            chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0/0/0\t0/1/1/1
            chr1\t300\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0/0/0\t0/0/0/0
            """
        ),
        encoding="utf-8",
    )
    popmap_path.write_text("s1\tPOP1\ns2\tPOP1\ns3\tPOP2\ns4\tPOP2\n", encoding="utf-8")

    loci = list(iter_loci_from_vcf(str(vcf_path), str(popmap_path)))
    rows, summary = calc_rho_windows(loci, window_size=1000, minimum_snps=1)

    assert len(rows) == 1
    assert rows[0]["pop1"] == "POP1"
    assert rows[0]["pop2"] == "POP2"
    assert summary["pairs"][0]["no_sites"] == 3


def test_iter_loci_from_vcf_skips_non_snp_records(tmp_path):
    vcf_path = tmp_path / "toy.vcf"
    popmap_path = tmp_path / "populations.tsv"

    vcf_path.write_text(
        textwrap.dedent(
            """\
            ##fileformat=VCFv4.2
            ##contig=<ID=chr1>
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
            chr1\t100\t.\tA\tAT\t.\tPASS\t.\tGT\t0/1\t0/0
            chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0
            chr1\t300\t.\tC\t.\t.\tPASS\t.\tGT\t0/0\t0/0
            """
        ),
        encoding="utf-8",
    )
    popmap_path.write_text("s1\tPOP1\ns2\tPOP1\n", encoding="utf-8")

    loci = list(iter_loci_from_vcf(str(vcf_path), str(popmap_path)))

    assert len(loci) == 2
    assert [record.position for record in [locus[0] for locus in loci]] == [200, 300]


def test_record_passes_ploidyscope_vcf_filters_matches_vcf_loader_rules():
    assert record_passes_ploidyscope_vcf_filters("A", ("G",)) is True
    assert record_passes_ploidyscope_vcf_filters("A", ()) is True
    assert record_passes_ploidyscope_vcf_filters("A", ("AT",)) is False
    assert record_passes_ploidyscope_vcf_filters("AT", ("A",)) is False
    assert record_passes_ploidyscope_vcf_filters("GTT", ("G", "*")) is False


def test_tajima_counts_fixed_alt_sites_in_watterson_theta():
    records = [
        ["POP1", "2", "chr1", "100", "2", "10", "0"],
        ["POP1", "2", "chr1", "200", "2", "10", "2"],
        ["POP1", "2", "chr1", "300", "4", "10", "0", "1"],
    ]

    tajima_rows = calc_tajima_d_windows(records, window_size=1000)

    assert len(tajima_rows) == 1
    assert tajima_rows[0]["raw_watterson_theta"] == pytest.approx(
        1.0 + (1.0 / (1.0 + 0.5 + (1.0 / 3.0)))
    )


def test_tajima_stdev_uses_pixy_mean_all_sites_convention():
    records = [
        ["POP1", "2", "chr1", "100", "4", "10", "0", "1"],
        ["POP1", "2", "chr1", "200", "4", "10", "1", "1"],
        ["POP1", "2", "chr1", "300", "4", "10", "-9", "-9"],
    ]

    rows = calc_tajima_d_windows(records, window_size=1000)

    assert len(rows) == 1
    assert rows[0]["no_sites"] == 2
    assert rows[0]["tajima_d_stdev"] == pytest.approx(0.0)
    assert rows[0]["tajima_d"] == "NA"


def test_tajima_stdev_is_na_when_window_has_no_observed_sites():
    records = [
        ["POP1", "2", "chr1", "100", "4", "10", "-9", "-9"],
        ["POP1", "2", "chr1", "200", "4", "10", "-9", "-9"],
    ]

    rows = calc_tajima_d_windows(records, window_size=1000)

    assert len(rows) == 1
    assert rows[0]["no_sites"] == 0
    assert rows[0]["raw_pi"] == 0.0
    assert rows[0]["raw_watterson_theta"] == 0.0
    assert rows[0]["tajima_d_stdev"] == "NA"
    assert rows[0]["tajima_d"] == "NA"


def test_tajima_stdev_is_na_when_mean_all_sites_rounds_below_two():
    records = [
        ["POP1", "2", "chr1", "100", "2", "10", "0"],
        ["POP1", "2", "chr1", "200", "2", "10", "-9"],
    ]

    rows = calc_tajima_d_windows(records, window_size=1000)

    assert len(rows) == 1
    assert rows[0]["no_sites"] == 1
    assert rows[0]["raw_pi"] == 0.0
    assert rows[0]["raw_watterson_theta"] == 0.0
    assert rows[0]["tajima_d_stdev"] == "NA"
    assert rows[0]["tajima_d"] == "NA"


def test_build_comparison_plan_routes_mixed_and_same_ploidy():
    samples = [
        ("a", "BL11-S2", 2, "non_serpentine"),
        ("b", "DKS", 2, "serpentine"),
        ("c", "TNA", 4, "non_serpentine"),
        ("d", "TDS", 4, "serpentine"),
    ]

    plans = build_comparison_plan(
        [
            SampleMetadata(sample=sample, population=pop, ploidy=ploidy, serpentine=serp)
            for sample, pop, ploidy, serp in samples
        ]
    )

    route_map = {(plan.pop1, plan.pop2): plan.route for plan in plans}
    assert route_map[("BL11-S2", "DKS")] == "pixy"
    assert route_map[("TDS", "TNA")] == "pixy"
    assert route_map[("BL11-S2", "TDS")] == "scantools"
    assert route_map[("DKS", "TNA")] == "scantools"


def test_parse_metadata_respects_selected_populations(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        textwrap.dedent(
            """\
            sample	prefix	sample_barcode	lane	ploidy	fq1	fq2	pop	serpentine
            bl11	p1	bl11	L1	2	f1	f2	BL11-S2	non_serpentine
            bl24	p2	bl24	L1	2	f1	f2	BL24-RED	serpentine
            tna1	p3	tna1	L1	4	f1	f2	TNA	non_serpentine
            """
        ),
        encoding="utf-8",
    )

    records = parse_metadata(metadata, {"BL11-S2", "TNA"})

    assert set(records) == {"bl11", "tna1"}
    assert records["bl11"].ploidy == 2
    assert records["tna1"].ploidy == 4


def test_parse_metadata_accepts_generic_population_header_and_optional_label(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        textwrap.dedent(
            """\
            sample	population	ploidy
            s1	POP1	2
            s2	POP2	4
            """
        ),
        encoding="utf-8",
    )

    records = parse_metadata(metadata)

    assert set(records) == {"s1", "s2"}
    assert records["s1"].population == "POP1"
    assert records["s2"].serpentine == "unknown"


def test_normalize_scantools_bpm_rows_aligns_zero_based_windows():
    rows = normalize_scantools_bpm_rows(
        [
            {
                "scaff": "chr1",
                "start": "0.0",
                "end": "10000.0",
                "Rho": "0.1",
                "Fst": "0.2",
                "dxy": "0.3",
            }
        ]
    )

    assert rows == [
        {
            "chromosome": "chr1",
            "window_pos_1": 1,
            "window_pos_2": 10000,
            "rho": "0.1",
            "fst": "0.2",
            "avg_dxy": "0.3",
        }
    ]


def test_write_svg_heatmap_creates_svg(tmp_path):
    output = tmp_path / "heatmap.svg"
    write_svg_heatmap(
        output,
        [
            {
                "comparison": "A__B",
                "route": "pixy",
                "metric": "dxy",
                "match": True,
                "rows_left": 10,
                "rows_right": 10,
                "shared_keys": 10,
                "mismatch_count": 0,
                "max_abs_diff": 0.0,
            },
            {
                "comparison": "A__B",
                "route": "pixy",
                "metric": "fst_wc",
                "match": False,
                "rows_left": 10,
                "rows_right": 8,
                "shared_keys": 8,
                "mismatch_count": 2,
                "max_abs_diff": 0.1,
            },
        ],
    )

    text = output.read_text(encoding="utf-8")
    assert text.startswith("<svg")
    assert "Serpentine Comparison Heatmap" in text
    assert "A__B" in text
