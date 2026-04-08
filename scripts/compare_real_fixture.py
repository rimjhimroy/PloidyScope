#!/usr/bin/env python3
"""Compare an existing fixture VCF against external baseline tools.

Usage:
    python scripts/compare_real_fixture.py

The helper assumes a fixture VCF plus metadata already exist and then:

1. Rebuilds helper files such as the population map and comparison map.
2. Routes same-ploidy comparisons to pixy and mixed-ploidy comparisons to
    ScanTools-style recoded tables.
3. Writes window-level comparison summaries against the original tools.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root))

from ploidyscope.stats.cli import write_tsv
from ploidyscope.stats.common import iter_loci_from_vcf
from ploidyscope.stats.diversity import calc_dxy_windows
from ploidyscope.stats.diversity import calc_dxy_windows_scantools_mean
from ploidyscope.stats.diversity import calc_hudson_fst_windows
from ploidyscope.stats.diversity import calc_pi_windows
from ploidyscope.stats.diversity import calc_tajima_d_windows
from ploidyscope.stats.rho import calc_rho_windows
from scripts.extract_fixture_region import DEFAULT_COMPARISON_DIR
from scripts.extract_fixture_region import DEFAULT_FIXTURE_DIR
from scripts.extract_fixture_region import DEFAULT_SITE_STRIDE
from scripts.extract_fixture_region import CANONICAL_OUTPUT_PREFIX
from scripts.extract_fixture_region import CANONICAL_REGION
from scripts.extract_fixture_region import PairPlan
from scripts.extract_fixture_region import SampleMetadata
from scripts.extract_fixture_region import build_comparison_plan
from scripts.extract_fixture_region import ensure_exists
from scripts.extract_fixture_region import index_vcfgz
from scripts.extract_fixture_region import open_vcf_writer
from scripts.extract_fixture_region import parse_metadata
from scripts.extract_fixture_region import read_vcf_samples
from scripts.extract_fixture_region import region_chromosome
from scripts.extract_fixture_region import write_minimal_vcf_header
from scripts.extract_fixture_region import write_minimal_vcf_record
from scripts.extract_fixture_region import write_comparison_map
from scripts.extract_fixture_region import write_lines
from scripts.extract_fixture_region import write_popmap


WC_FST_NOTE = "Weir-Cockerham fst is only compared for diploid populations because pixy 2.0.0.beta14 restricts fst to diploid genotypes. For fixture validation, diploid WC uses the same filtered site class as PloidyScope's VCF loader and compares only non-NA windows; the internal diploid WC accumulation now follows pixy's separate-nansum window semantics."


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metadata",
        default=str(DEFAULT_FIXTURE_DIR / f"{CANONICAL_OUTPUT_PREFIX}_metadata.tsv"),
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_COMPARISON_DIR),
    )
    parser.add_argument(
        "--fixture-vcf",
        default=str(DEFAULT_FIXTURE_DIR / f"{CANONICAL_OUTPUT_PREFIX}.vcf.gz"),
    )
    parser.add_argument(
        "--fixture-popmap",
        default=str(DEFAULT_FIXTURE_DIR / f"{CANONICAL_OUTPUT_PREFIX}_popmap.tsv"),
    )
    parser.add_argument(
        "--comparison-map",
        default=str(DEFAULT_FIXTURE_DIR / "comparison_map.tsv"),
    )
    parser.add_argument(
        "--region",
        default=CANONICAL_REGION,
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "--site-stride",
        type=int,
        default=DEFAULT_SITE_STRIDE,
    )
    parser.add_argument(
        "--minimum-snps",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--pixy-python",
        default=str(repo_root / ".tool-baselines" / "pixy-venv" / "bin" / "python"),
    )
    parser.add_argument(
        "--pixy-exe",
        default=str(repo_root / ".tool-baselines" / "pixy-venv" / "bin" / "pixy"),
    )
    parser.add_argument(
        "--pixy-cores",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--tabix-exe",
        default=str(repo_root / ".tool-baselines" / "pixy-venv" / "bin" / "tabix"),
    )
    parser.add_argument(
        "--scantools-bpm",
        default=str(repo_root / ".tool-baselines" / "scantools" / "bpm.py"),
    )
    return parser


def run(command: list[str], env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env=env,
    )


def write_subset_vcf(input_vcf: Path, output_vcf: Path, sample_names: list[str]) -> None:
    try:
        from pysam import VariantFile
    except ModuleNotFoundError:
        VariantFile = None

    if VariantFile is not None and output_vcf.suffix == ".gz":
        with VariantFile(str(input_vcf)) as source:
            source.subset_samples(sample_names)
            contig_names = list(source.header.contigs)
            primary_contig = contig_names[0] if contig_names else "unknown"
            contig_length = source.header.contigs[primary_contig].length if contig_names else None
            with open_vcf_writer(output_vcf) as target:
                write_minimal_vcf_header(target, primary_contig, contig_length, sample_names)
                for record in source:
                    write_minimal_vcf_record(target, record, sample_names)
        return

    keep: list[int] | None = None
    open_input = gzip.open if input_vcf.suffix == ".gz" else open
    with open_input(input_vcf, "rt", encoding="utf-8") as source, output_vcf.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith("##"):
                target.write(line)
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                header_samples = header[9:]
                keep = [index for index, sample in enumerate(header_samples) if sample in sample_names]
                target.write("\t".join(header[:9] + [header_samples[index] for index in keep]) + "\n")
                continue
            fields = line.rstrip("\n").split("\t")
            assert keep is not None
            target.write("\t".join(fields[:9] + [fields[9 + index] for index in keep]) + "\n")


def record_passes_ploidyscope_vcf_filters(ref_allele: str, alternates: tuple[str, ...] | None) -> bool:
    alternates = alternates or ()
    if len(alternates) > 1:
        return False
    if len(str(ref_allele)) != 1:
        return False
    if alternates and any(len(str(allele)) != 1 for allele in alternates):
        return False
    return True


def write_ploidyscope_filtered_vcf(input_vcf: Path, output_vcf: Path) -> None:
    try:
        from pysam import VariantFile
    except ModuleNotFoundError:
        VariantFile = None

    if VariantFile is not None:
        with VariantFile(str(input_vcf)) as source:
            sample_names = list(source.header.samples)
            contig_names = list(source.header.contigs)
            primary_contig = contig_names[0] if contig_names else "unknown"
            contig_length = source.header.contigs[primary_contig].length if contig_names else None
            with open_vcf_writer(output_vcf) as target:
                write_minimal_vcf_header(target, primary_contig, contig_length, sample_names)
                for record in source:
                    if not record_passes_ploidyscope_vcf_filters(record.ref, record.alts):
                        continue
                    write_minimal_vcf_record(target, record, sample_names)
        return

    open_input = gzip.open if input_vcf.suffix == ".gz" else open
    with open_input(input_vcf, "rt", encoding="utf-8") as source, open_vcf_writer(output_vcf) as target:
        for line in source:
            if line.startswith("##") or line.startswith("#CHROM"):
                target.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            ref_allele = fields[3]
            alternates = () if fields[4] == "." else tuple(fields[4].split(","))
            if not record_passes_ploidyscope_vcf_filters(ref_allele, alternates):
                continue
            target.write(line)


def build_strict_scantools_table_from_loci(loci: Iterable[list[object]]) -> list[list[str]]:
    strict_rows: list[list[str]] = []
    for locus in loci:
        if all(sum(genotype is not None for genotype in record.genotypes) >= 2 for record in locus):
            for record in locus:
                row = [
                    record.population,
                    str(record.ploidy),
                    record.chromosome,
                    str(record.position),
                    str(record.possible_alleles),
                    "0",
                ]
                row.extend("-9" if genotype is None else str(genotype) for genotype in record.genotypes)
                strict_rows.append(row)
    return strict_rows


def first_data_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: object) -> float | None:
    text = str(value)
    if text in {"NA", "-9", "unsupported", "None"}:
        return None
    return float(text)


def compare_tables(
    left_rows: list[dict[str, object]],
    right_rows: list[dict[str, object]],
    left_key_fields: list[str],
    right_key_fields: list[str],
    left_value_field: str,
    right_value_field: str,
    abs_tol: float = 1e-10,
) -> dict[str, object]:
    left_index = {tuple(str(row[field]) for field in left_key_fields): row for row in left_rows}
    right_index = {tuple(str(row[field]) for field in right_key_fields): row for row in right_rows}
    shared_keys = sorted(set(left_index) & set(right_index))
    only_left = sorted(set(left_index) - set(right_index))
    only_right = sorted(set(right_index) - set(left_index))

    mismatch_examples: list[dict[str, object]] = []
    mismatch_count = 0
    max_abs_diff = 0.0
    compared_rows = 0
    for key in shared_keys:
        left_value = as_float(left_index[key][left_value_field])
        right_value = as_float(right_index[key][right_value_field])
        if left_value is None or right_value is None:
            if str(left_index[key][left_value_field]) != str(right_index[key][right_value_field]):
                mismatch_count += 1
                if len(mismatch_examples) < 5:
                    mismatch_examples.append(
                        {
                            "key": key,
                            "left": left_index[key][left_value_field],
                            "right": right_index[key][right_value_field],
                        }
                    )
            compared_rows += 1
            continue
        difference = abs(left_value - right_value)
        max_abs_diff = max(max_abs_diff, difference)
        compared_rows += 1
        if difference > abs_tol:
            mismatch_count += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append(
                    {
                        "key": key,
                        "left": left_value,
                        "right": right_value,
                        "difference": left_value - right_value,
                    }
                )

    return {
        "rows_left": len(left_rows),
        "rows_right": len(right_rows),
        "shared_keys": len(shared_keys),
        "only_left": len(only_left),
        "only_right": len(only_right),
        "compared_rows": compared_rows,
        "mismatch_count": mismatch_count,
        "max_abs_diff": max_abs_diff,
        "match": mismatch_count == 0 and not only_left and not only_right,
        "abs_tol": abs_tol,
        "mismatch_examples": mismatch_examples,
    }


def normalize_scantools_bpm_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for row in rows:
        if row["scaff"] == "Genome":
            continue
        normalized.append(
            {
                "chromosome": row["scaff"],
                "window_pos_1": int(float(row["start"])) + 1,
                "window_pos_2": int(float(row["end"])),
                "rho": row["Rho"],
                "fst": row["Fst"],
                "avg_dxy": row["dxy"],
            }
        )
    return normalized


def comparison_status(result: dict[str, object]) -> str:
    match = result.get("match")
    if match is True:
        return "exact"
    if match == "unsupported":
        return "unsupported"
    if match == "baseline_only":
        return "baseline_only"
    if int(result.get("shared_keys", 0)) == 0:
        return "no_shared_windows"
    if int(result.get("mismatch_count", 0)) == 0:
        if int(result.get("only_left", 0)) == 0 and int(result.get("only_right", 0)) == 0:
            return "exact"
        return "key_mismatch_only"
    return "value_diff"


def coverage_bar(shared_keys: int, rows_left: int, rows_right: int, width: int = 12) -> str:
    total = max(rows_left, rows_right, 1)
    filled = int(round((shared_keys / total) * width))
    filled = max(0, min(width, filled))
    return ("#" * filled) + ("." * (width - filled))


def status_fill(status: str) -> str:
    fills = {
        "exact": "#2e8b57",
        "key_mismatch_only": "#d4a72c",
        "value_diff": "#c95a49",
        "unsupported": "#7f8c8d",
        "baseline_only": "#4c78a8",
        "no_shared_windows": "#9c755f",
    }
    return fills.get(status, "#b0b0b0")


def metric_sort_key(metric: str) -> tuple[int, str]:
    metric_order = {
        "pi": 0,
        "raw_pi": 1,
        "watterson_theta": 2,
        "tajima_d": 3,
        "tajima_d_stdev": 4,
        "dxy": 5,
        "dxy_weighted": 6,
        "dxy_scantools_mean": 7,
        "fst_wc": 8,
        "fst_hudson": 9,
        "fst": 10,
        "rho": 11,
    }
    return metric_order.get(metric, 99), metric


def write_svg_heatmap(path: Path, summary_rows: list[dict[str, object]]) -> None:
    ordered_rows = sorted(summary_rows, key=lambda row: (str(row["route"]), str(row["comparison"]), metric_sort_key(str(row["metric"]))))
    comparisons = []
    for row in ordered_rows:
        comparison = str(row["comparison"])
        if comparison not in comparisons:
            comparisons.append(comparison)
    metrics = sorted({str(row["metric"]) for row in ordered_rows}, key=metric_sort_key)
    row_lookup = {(str(row["comparison"]), str(row["metric"])): row for row in ordered_rows}

    cell_width = 90
    cell_height = 26
    left_margin = 180
    top_margin = 80
    legend_height = 72
    width = left_margin + (len(metrics) * cell_width) + 40
    height = top_margin + (len(comparisons) * cell_height) + legend_height + 40

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfbf8"/>',
        '<text x="24" y="30" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2933">Fixture Comparison Heatmap</text>',
        '<text x="24" y="52" font-family="Arial, sans-serif" font-size="11" fill="#52606d">Green = exact, red = value difference, yellow = key mismatch, gray = unsupported.</text>',
    ]

    for column_index, metric in enumerate(metrics):
        x = left_margin + (column_index * cell_width) + (cell_width / 2)
        parts.append(
            f'<text x="{x}" y="68" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#1f2933">{html.escape(metric)}</text>'
        )

    for row_index, comparison in enumerate(comparisons):
        y = top_margin + (row_index * cell_height)
        parts.append(
            f'<text x="{left_margin - 8}" y="{y + 17}" text-anchor="end" font-family="Arial, sans-serif" font-size="10" fill="#1f2933">{html.escape(comparison)}</text>'
        )
        for column_index, metric in enumerate(metrics):
            x = left_margin + (column_index * cell_width)
            row = row_lookup.get((comparison, metric))
            if row is None:
                fill = "#ececec"
                label = "-"
                title = f"{comparison} | {metric} | missing"
            else:
                status = comparison_status(row)
                fill = status_fill(status)
                label = {
                    "exact": "E",
                    "key_mismatch_only": "K",
                    "value_diff": "V",
                    "unsupported": "U",
                    "baseline_only": "B",
                    "no_shared_windows": "N",
                }.get(status, "?")
                title = (
                    f"{comparison} | {metric} | status={status} | shared={row['shared_keys']} | "
                    f"mismatches={row['mismatch_count']} | max_abs_diff={row['max_abs_diff']}"
                )
            parts.append(
                f'<g><title>{html.escape(title)}</title><rect x="{x}" y="{y}" width="{cell_width - 4}" height="{cell_height - 4}" rx="4" ry="4" fill="{fill}" stroke="#ffffff" stroke-width="1"/><text x="{x + (cell_width - 4) / 2}" y="{y + 16}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" font-weight="700" fill="#ffffff">{label}</text></g>'
            )

    legend_items = [
        ("exact", "E"),
        ("key_mismatch_only", "K"),
        ("value_diff", "V"),
        ("unsupported", "U"),
        ("baseline_only", "B"),
        ("no_shared_windows", "N"),
    ]
    legend_y = top_margin + (len(comparisons) * cell_height) + 24
    for index, (status, label) in enumerate(legend_items):
        x = 24 + (index * 125)
        fill = status_fill(status)
        parts.append(f'<rect x="{x}" y="{legend_y}" width="18" height="18" rx="3" ry="3" fill="{fill}"/>')
        parts.append(f'<text x="{x + 9}" y="{legend_y + 13}" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" font-weight="700" fill="#ffffff">{label}</text>')
        parts.append(f'<text x="{x + 26}" y="{legend_y + 13}" font-family="Arial, sans-serif" font-size="10" fill="#1f2933">{html.escape(status)}</text>')

    parts.append('</svg>')
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def write_visual_summary(
    path: Path,
    summary_rows: list[dict[str, object]],
    comparison_plan: list[PairPlan],
    selected_pops: list[str],
    excluded_pops: list[str],
    region: str,
    window_size: int,
    site_stride: int,
) -> None:
    pixy_pairs = [f"{plan.pop1}__{plan.pop2}" for plan in comparison_plan if plan.route == "pixy"]
    scantools_pairs = [f"{plan.pop1}__{plan.pop2}" for plan in comparison_plan if plan.route == "scantools"]
    ordered_rows = sorted(summary_rows, key=lambda row: (str(row["route"]), str(row["comparison"]), str(row["metric"])))

    lines = [
        "# Fixture Comparison Summary",
        "",
        "![SVG heatmap](visual_summary.svg)",
        "",
        f"Region: `{region}`  ",
        f"Window size: `{window_size}`  ",
        f"Site stride: `{site_stride}`  ",
        f"Selected populations: `{', '.join(selected_pops)}`  ",
        f"Excluded populations: `{', '.join(excluded_pops) if excluded_pops else 'none'}`",
        "",
        "```mermaid",
        "flowchart TD",
        "    Fixture[Serpentine 10 Mb Fixture]",
        f"    Fixture --> Diploid[Diploid pops: {' / '.join(sorted(pop for pop in selected_pops if pop in {'BL11-S2', 'DKS'}))}]",
        f"    Fixture --> Tetraploid[Tetraploid pops: {' / '.join(sorted(pop for pop in selected_pops if pop in {'TDS', 'TNA'}))}]",
        f"    Diploid --> PixyDip[Pixy route: {'; '.join(pixy_pairs[:1]) or 'none'}]",
        f"    Tetraploid --> PixyTet[Pixy route: {'; '.join(pixy_pairs[1:]) or 'none'}]",
        f"    Fixture --> Mixed[ScanTools route: {'; '.join(scantools_pairs) or 'none'}]",
        "```",
        "",
        "| comparison | route | metric | status | shared | coverage | max abs diff |",
        "| --- | --- | --- | --- | ---: | --- | ---: |",
    ]
    for row in ordered_rows:
        shared = int(row["shared_keys"]) if str(row["shared_keys"]) != "NA" else 0
        rows_left = int(row["rows_left"]) if str(row["rows_left"]) != "NA" else 0
        rows_right = int(row["rows_right"]) if str(row["rows_right"]) != "NA" else 0
        diff = row["max_abs_diff"]
        lines.append(
            "| {comparison} | {route} | {metric} | {status} | {shared} | {bar} | {diff} |".format(
                comparison=row["comparison"],
                route=row["route"],
                metric=row["metric"],
                status=comparison_status(row),
                shared=shared,
                bar=coverage_bar(shared, rows_left, rows_right),
                diff=diff,
            )
        )

    write_lines(path, lines)


def rho_rows_to_fst_rows(rho_rows: list[dict[str, object]], minimum_snps: int) -> list[dict[str, object]]:
    fst_rows: list[dict[str, object]] = []
    for row in rho_rows:
        value: object = "NA"
        fst_no_snps = int(row.get("fst_no_snps", row["no_snps"]))
        if fst_no_snps >= minimum_snps and float(row["fst_den"]) != 0.0:
            value = float(row["fst_num"]) / float(row["fst_den"])
        fst_rows.append(
            {
                "pop1": row["pop1"],
                "pop2": row["pop2"],
                "chromosome": row["chromosome"],
                "window_pos_1": row["window_pos_1"],
                "window_pos_2": row["window_pos_2"],
                "fst": value,
            }
        )
    return fst_rows


def group_rows(rows: list[dict[str, object]], field: str, value: str) -> list[dict[str, object]]:
    return [row for row in rows if str(row[field]) == value]


def summarize_to_tsv(summary_rows: list[dict[str, object]], path: Path) -> None:
    write_tsv(summary_rows, str(path))


def run_pixy_stats(
    pixy_exe: Path,
    pixy_env: dict[str, str],
    split_vcfgz: Path,
    popmap: Path,
    output_dir: Path,
    output_prefix: str,
    stats: list[str],
    window_size: int,
    n_cores: int,
    fst_type: str | None = None,
) -> None:
    command = [
        str(pixy_exe),
        "--stats",
        *stats,
        "--vcf",
        str(split_vcfgz),
        "--populations",
        str(popmap),
        "--window_size",
        str(window_size),
        "--output_folder",
        str(output_dir),
        "--output_prefix",
        output_prefix,
        "--n_cores",
        str(n_cores),
        "--silent",
    ]
    if fst_type is not None:
        command.extend(["--fst_type", fst_type])
    run(command, env=pixy_env)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    fixture_vcfgz = Path(args.fixture_vcf)
    metadata_path = Path(args.metadata)
    fixture_popmap = Path(args.fixture_popmap)
    comparison_map = Path(args.comparison_map)
    outdir = Path(args.outdir)
    pixy_python = Path(args.pixy_python)
    pixy_exe = Path(args.pixy_exe)
    tabix_exe = Path(args.tabix_exe)
    scantools_bpm = Path(args.scantools_bpm)

    ensure_exists(fixture_vcfgz, "fixture VCF")
    ensure_exists(metadata_path, "fixture metadata")
    ensure_exists(pixy_python, "pixy python")
    ensure_exists(pixy_exe, "pixy executable")
    ensure_exists(tabix_exe, "tabix executable")
    ensure_exists(scantools_bpm, "ScanTools BPM script")

    metadata_records = parse_metadata(metadata_path)
    vcf_samples = set(read_vcf_samples(fixture_vcfgz))
    selected_samples = sorted(
        sample for sample in metadata_records if sample in vcf_samples
    )
    missing_samples = sorted(sample for sample in metadata_records if sample not in vcf_samples)
    if missing_samples:
        raise ValueError(f"Selected metadata samples missing from the VCF: {', '.join(missing_samples[:10])}")
    selected_metadata = [metadata_records[sample] for sample in selected_samples]
    comparison_plan = build_comparison_plan(selected_metadata)
    selected_pops = sorted({sample.population for sample in selected_metadata})

    fixture_dir = outdir / "fixture"
    pixy_dir = outdir / "pixy"
    scantools_dir = outdir / "scantools"
    ploidyscope_dir = outdir / "ploidyscope"
    fixture_dir.mkdir(parents=True, exist_ok=True)
    for generated_dir in (pixy_dir, scantools_dir, ploidyscope_dir):
        if generated_dir.exists():
            shutil.rmtree(generated_dir)
        generated_dir.mkdir(parents=True, exist_ok=True)

    write_popmap(fixture_popmap, selected_metadata)
    write_comparison_map(comparison_map, comparison_plan)

    samples_by_population: dict[str, list[SampleMetadata]] = defaultdict(list)
    samples_by_ploidy: dict[int, list[SampleMetadata]] = defaultdict(list)
    for sample in selected_metadata:
        samples_by_population[sample.population].append(sample)
        samples_by_ploidy[sample.ploidy].append(sample)

    pixy_env = os.environ.copy()
    pixy_env["PATH"] = str(pixy_exe.parent) + os.pathsep + pixy_env.get("PATH", "")

    summary_rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "fixture": {
            "vcf": str(fixture_vcfgz),
            "metadata": str(metadata_path),
            "popmap": str(fixture_popmap),
            "comparison_map": str(comparison_map),
            "region": args.region,
            "chromosome": region_chromosome(args.region),
            "window_size": args.window_size,
            "site_stride": args.site_stride,
            "minimum_snps": args.minimum_snps,
            "selected_populations": selected_pops,
            "excluded_populations": [],
            "notes": [
                WC_FST_NOTE,
                "The comparison helper expects a pre-extracted fixture VCF and rewrites only helper files plus regenerated baseline outputs.",
            ],
        },
        "plans": [asdict(plan) for plan in comparison_plan],
        "pixy": {},
        "scantools": {},
    }

    for ploidy in sorted(samples_by_ploidy):
        samples = sorted(samples_by_ploidy[ploidy], key=lambda item: (item.population, item.sample))
        populations = sorted({sample.population for sample in samples})
        split_dir = pixy_dir / f"ploidy_{ploidy}"
        split_dir.mkdir(parents=True, exist_ok=True)

        split_vcfgz = split_dir / f"serpentine_10mb.ploidy_{ploidy}.vcf.gz"
        split_popmap = split_dir / f"serpentine_10mb.ploidy_{ploidy}.popmap.tsv"
        write_subset_vcf(fixture_vcfgz, split_vcfgz, [sample.sample for sample in samples])
        index_vcfgz(split_vcfgz, pixy_python)
        write_popmap(split_popmap, samples)

        prefix = f"serpentine_10mb_ploidy_{ploidy}"
        run_pixy_stats(
            pixy_exe=pixy_exe,
            pixy_env=pixy_env,
            split_vcfgz=split_vcfgz,
            popmap=split_popmap,
            output_dir=split_dir,
            output_prefix=prefix,
            stats=["pi", "watterson_theta", "tajima_d", "dxy"],
            window_size=args.window_size,
            n_cores=args.pixy_cores,
        )
        run_pixy_stats(
            pixy_exe=pixy_exe,
            pixy_env=pixy_env,
            split_vcfgz=split_vcfgz,
            popmap=split_popmap,
            output_dir=split_dir,
            output_prefix=f"{prefix}.hudson",
            stats=["fst"],
            window_size=args.window_size,
            n_cores=args.pixy_cores,
            fst_type="hudson",
        )
        if ploidy == 2:
            split_wc_vcfgz = split_dir / f"serpentine_10mb.ploidy_{ploidy}.wc_compatible.vcf.gz"
            write_ploidyscope_filtered_vcf(split_vcfgz, split_wc_vcfgz)
            index_vcfgz(split_wc_vcfgz, pixy_python)
            run_pixy_stats(
                pixy_exe=pixy_exe,
                pixy_env=pixy_env,
                split_vcfgz=split_wc_vcfgz,
                popmap=split_popmap,
                output_dir=split_dir,
                output_prefix=prefix,
                stats=["fst"],
                window_size=args.window_size,
                n_cores=args.pixy_cores,
                fst_type="wc",
            )

        loci = list(iter_loci_from_vcf(str(split_vcfgz), str(split_popmap)))
        pi_rows = calc_pi_windows(loci, window_size=args.window_size)
        tajima_rows = calc_tajima_d_windows(loci, window_size=args.window_size)
        dxy_rows = calc_dxy_windows(loci, window_size=args.window_size)
        hudson_rows = calc_hudson_fst_windows(
            loci,
            window_size=args.window_size,
            minimum_snps=args.minimum_snps,
        )
        rho_rows, _ = calc_rho_windows(loci, window_size=args.window_size, minimum_snps=args.minimum_snps)
        fst_wc_rows = rho_rows_to_fst_rows(rho_rows, args.minimum_snps)

        write_tsv(pi_rows, str(ploidyscope_dir / f"ploidy_{ploidy}_pi.tsv"))
        write_tsv(tajima_rows, str(ploidyscope_dir / f"ploidy_{ploidy}_tajima_d.tsv"))
        write_tsv(dxy_rows, str(ploidyscope_dir / f"ploidy_{ploidy}_dxy.tsv"))
        write_tsv(fst_wc_rows, str(ploidyscope_dir / f"ploidy_{ploidy}_fst_wc.tsv"))
        write_tsv(hudson_rows, str(ploidyscope_dir / f"ploidy_{ploidy}_fst_hudson.tsv"))

        pixy_pi_rows = first_data_rows(split_dir / f"{prefix}_pi.txt")
        pixy_theta_rows = first_data_rows(split_dir / f"{prefix}_watterson_theta.txt")
        pixy_tajima_rows = first_data_rows(split_dir / f"{prefix}_tajima_d.txt")
        pixy_dxy_rows = first_data_rows(split_dir / f"{prefix}_dxy.txt")

        ploidy_summary: dict[str, object] = {
            "populations": populations,
            "population_stats": {},
        }
        for population in populations:
            pi_compare = compare_tables(
                left_rows=group_rows(pi_rows, "pop", population),
                right_rows=group_rows(pixy_pi_rows, "pop", population),
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="pi",
                right_value_field="avg_pi",
            )
            theta_compare = compare_tables(
                left_rows=group_rows(tajima_rows, "pop", population),
                right_rows=group_rows(pixy_theta_rows, "pop", population),
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="raw_watterson_theta",
                right_value_field="raw_watterson_theta",
            )
            tajima_compare = compare_tables(
                left_rows=group_rows(tajima_rows, "pop", population),
                right_rows=group_rows(pixy_tajima_rows, "pop", population),
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="tajima_d",
                right_value_field="tajima_d",
            )
            raw_pi_compare = compare_tables(
                left_rows=group_rows(tajima_rows, "pop", population),
                right_rows=group_rows(pixy_tajima_rows, "pop", population),
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="raw_pi",
                right_value_field="raw_pi",
            )
            stdev_compare = compare_tables(
                left_rows=group_rows(tajima_rows, "pop", population),
                right_rows=group_rows(pixy_tajima_rows, "pop", population),
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="tajima_d_stdev",
                right_value_field="tajima_d_stdev",
            )
            ploidy_summary["population_stats"][population] = {
                "pi": pi_compare,
                "watterson_theta": theta_compare,
                "tajima_d": tajima_compare,
                "raw_pi": raw_pi_compare,
                "tajima_d_stdev": stdev_compare,
            }
            for metric, result in ploidy_summary["population_stats"][population].items():
                summary_rows.append(
                    {
                        "comparison": population,
                        "route": "pixy",
                        "metric": metric,
                        "match": result["match"],
                        "rows_left": result["rows_left"],
                        "rows_right": result["rows_right"],
                        "shared_keys": result["shared_keys"],
                        "mismatch_count": result["mismatch_count"],
                        "max_abs_diff": result["max_abs_diff"],
                        "note": "within-population same-ploidy comparison",
                    }
                )

        same_ploidy_plans = [
            plan for plan in comparison_plan if plan.route == "pixy" and plan.ploidy1 == ploidy and plan.ploidy2 == ploidy
        ]
        if same_ploidy_plans:
            plan = same_ploidy_plans[0]
            pair_name = f"{plan.pop1}__{plan.pop2}"
            pair_dxy_compare = compare_tables(
                left_rows=dxy_rows,
                right_rows=pixy_dxy_rows,
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="avg_dxy",
                right_value_field="avg_dxy",
            )
            pixy_hudson_rows = first_data_rows(split_dir / f"{prefix}.hudson_fst.txt")
            fst_hudson_compare = compare_tables(
                left_rows=[row for row in hudson_rows if row["avg_hudson_fst"] != "NA"],
                right_rows=pixy_hudson_rows,
                left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                left_value_field="avg_hudson_fst",
                right_value_field="avg_hudson_fst",
            )
            ploidy_summary[pair_name] = {
                "dxy": pair_dxy_compare,
                "fst_hudson": fst_hudson_compare,
            }
            summary_rows.append(
                {
                    "comparison": pair_name,
                    "route": "pixy",
                    "metric": "dxy",
                    "match": pair_dxy_compare["match"],
                    "rows_left": pair_dxy_compare["rows_left"],
                    "rows_right": pair_dxy_compare["rows_right"],
                    "shared_keys": pair_dxy_compare["shared_keys"],
                    "mismatch_count": pair_dxy_compare["mismatch_count"],
                    "max_abs_diff": pair_dxy_compare["max_abs_diff"],
                    "note": "same-ploidy pairwise dxy",
                }
            )
            summary_rows.append(
                {
                    "comparison": pair_name,
                    "route": "pixy",
                    "metric": "fst_hudson",
                    "match": fst_hudson_compare["match"],
                    "rows_left": fst_hudson_compare["rows_left"],
                    "rows_right": fst_hudson_compare["rows_right"],
                    "shared_keys": fst_hudson_compare["shared_keys"],
                    "mismatch_count": fst_hudson_compare["mismatch_count"],
                    "max_abs_diff": fst_hudson_compare["max_abs_diff"],
                    "note": "Hudson fst compared against pixy on same-ploidy split VCFs",
                }
            )
            if ploidy == 2:
                pixy_fst_wc_rows = [
                    row for row in first_data_rows(split_dir / f"{prefix}_fst.txt") if row["avg_wc_fst"] != "NA"
                ]
                fst_wc_compare = compare_tables(
                    left_rows=[row for row in fst_wc_rows if row["fst"] != "NA"],
                    right_rows=pixy_fst_wc_rows,
                    left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                    right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
                    left_value_field="fst",
                    right_value_field="avg_wc_fst",
                )
                ploidy_summary[pair_name]["fst_wc"] = fst_wc_compare
                summary_rows.append(
                    {
                        "comparison": pair_name,
                        "route": "pixy",
                        "metric": "fst_wc",
                        "match": fst_wc_compare["match"],
                        "rows_left": fst_wc_compare["rows_left"],
                        "rows_right": fst_wc_compare["rows_right"],
                        "shared_keys": fst_wc_compare["shared_keys"],
                        "mismatch_count": fst_wc_compare["mismatch_count"],
                        "max_abs_diff": fst_wc_compare["max_abs_diff"],
                        "note": WC_FST_NOTE,
                    }
                )
            else:
                ploidy_summary[pair_name]["fst_hudson"] = fst_hudson_compare
                summary_rows.append(
                    {
                        "comparison": pair_name,
                        "route": "pixy",
                        "metric": "fst_wc",
                        "match": "unsupported",
                        "rows_left": 0,
                        "rows_right": 0,
                        "shared_keys": 0,
                        "mismatch_count": 0,
                        "max_abs_diff": "NA",
                        "note": WC_FST_NOTE,
                    }
                )

        summary["pixy"][str(ploidy)] = ploidy_summary

    for plan in comparison_plan:
        if plan.route != "scantools":
            continue
        pair_name = f"{plan.pop1}__{plan.pop2}"
        pair_dir = scantools_dir / pair_name
        pair_dir.mkdir(parents=True, exist_ok=True)
        pair_loci = list(
            iter_loci_from_vcf(
                str(fixture_vcfgz),
                str(fixture_popmap),
                populations=[plan.pop1, plan.pop2],
            )
        )
        strict_rows = build_strict_scantools_table_from_loci(pair_loci)
        strict_path = pair_dir / f"{pair_name}.strict.scantools.tsv"
        write_lines(strict_path, ["\t".join(row) for row in strict_rows])

        rho_rows, _ = calc_rho_windows(strict_rows, window_size=args.window_size, minimum_snps=args.minimum_snps)
        dxy_rows = calc_dxy_windows(strict_rows, window_size=args.window_size)
        dxy_scantools_mean_rows = calc_dxy_windows_scantools_mean(strict_rows, window_size=args.window_size)
        fst_rows = rho_rows_to_fst_rows(rho_rows, args.minimum_snps)
        write_tsv(rho_rows, str(ploidyscope_dir / f"{pair_name}_rho.tsv"))
        write_tsv(dxy_rows, str(ploidyscope_dir / f"{pair_name}_dxy.tsv"))
        write_tsv(dxy_scantools_mean_rows, str(ploidyscope_dir / f"{pair_name}_dxy_scantools_mean.tsv"))
        write_tsv(fst_rows, str(ploidyscope_dir / f"{pair_name}_fst.tsv"))

        scan_prefix = pair_name
        run(
            [
                str(pixy_python),
                str(scantools_bpm),
                "-i",
                str(strict_path),
                "-o",
                str(pair_dir) + os.sep,
                "-prefix",
                scan_prefix,
                "-ws",
                str(args.window_size),
                "-ms",
                str(args.minimum_snps),
                "-np",
                "2",
            ]
        )
        scan_rows = normalize_scantools_bpm_rows(first_data_rows(pair_dir / f"{scan_prefix}_BPM.txt"))

        rho_compare = compare_tables(
            left_rows=rho_rows,
            right_rows=scan_rows,
            left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            left_value_field="rho",
            right_value_field="rho",
        )
        fst_compare = compare_tables(
            left_rows=fst_rows,
            right_rows=scan_rows,
            left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            left_value_field="fst",
            right_value_field="fst",
        )
        dxy_weighted_compare = compare_tables(
            left_rows=dxy_rows,
            right_rows=scan_rows,
            left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            left_value_field="avg_dxy",
            right_value_field="avg_dxy",
        )
        dxy_scantools_mean_compare = compare_tables(
            left_rows=dxy_scantools_mean_rows,
            right_rows=scan_rows,
            left_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            right_key_fields=["chromosome", "window_pos_1", "window_pos_2"],
            left_value_field="dxy_scantools_mean",
            right_value_field="avg_dxy",
        )
        summary["scantools"][pair_name] = {
            "rho": rho_compare,
            "fst": fst_compare,
            "dxy_weighted": dxy_weighted_compare,
            "dxy_scantools_mean": dxy_scantools_mean_compare,
        }
        for metric, result in summary["scantools"][pair_name].items():
            summary_rows.append(
                {
                    "comparison": pair_name,
                    "route": "scantools",
                    "metric": metric,
                    "match": result["match"],
                    "rows_left": result["rows_left"],
                    "rows_right": result["rows_right"],
                    "shared_keys": result["shared_keys"],
                    "mismatch_count": result["mismatch_count"],
                    "max_abs_diff": result["max_abs_diff"],
                    "note": (
                        "mixed-ploidy pairwise comparison routed through strict recoding"
                        if metric in {"rho", "fst"}
                        else (
                            "mixed-ploidy weighted dxy compared against ScanTools site-mean dxy"
                            if metric == "dxy_weighted"
                            else "mixed-ploidy ScanTools-style site-mean dxy compared against ScanTools"
                        )
                    ),
                }
            )

    summary_json = outdir / "summary.json"
    summary_tsv = outdir / "summary.tsv"
    visual_summary = outdir / "visual_summary.md"
    visual_svg = outdir / "visual_summary.svg"
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summarize_to_tsv(summary_rows, summary_tsv)
    write_visual_summary(
        path=visual_summary,
        summary_rows=summary_rows,
        comparison_plan=comparison_plan,
        selected_pops=selected_pops,
        excluded_pops=[],
        region=args.region,
        window_size=args.window_size,
        site_stride=args.site_stride,
    )
    write_svg_heatmap(visual_svg, summary_rows)

    print(f"Summary written to {summary_json}")
    print(f"Tabular summary written to {summary_tsv}")
    print(f"Visual summary written to {visual_summary}")
    print(f"SVG heatmap written to {visual_svg}")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()