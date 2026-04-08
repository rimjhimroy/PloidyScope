from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from math import sqrt
from typing import Iterable, Optional, Sequence

from .common import PopulationRecord
from .common import harmonic_number
from .common import iter_loci
from .common import second_harmonic_number
from .common import window_bounds


def _choose_two(value: int) -> int:
    if value < 2:
        return 0
    return value * (value - 1) // 2


def _site_pi(record: PopulationRecord) -> tuple[int, int, int, float, bool]:
    observed = record.observed_alleles
    possible = record.possible_alleles
    alt = record.observed_alt
    ref = observed - alt
    diffs = alt * ref
    comparisons = _choose_two(observed)
    missing = _choose_two(possible) - comparisons
    raw_pi = (diffs / comparisons) if comparisons else 0.0
    is_variant = alt > 0 and ref > 0
    return diffs, comparisons, missing, raw_pi, is_variant


def _site_dxy(record1: PopulationRecord, record2: PopulationRecord) -> tuple[int, int, int]:
    ref1 = record1.observed_alleles - record1.observed_alt
    ref2 = record2.observed_alleles - record2.observed_alt
    diffs = (ref1 * record2.observed_alt) + (record1.observed_alt * ref2)
    comparisons = record1.observed_alleles * record2.observed_alleles
    missing = (record1.possible_alleles * record2.possible_alleles) - comparisons
    return diffs, comparisons, missing


def _site_hudson_fst(
    record1: PopulationRecord,
    record2: PopulationRecord,
) -> tuple[float, float] | None:
    allele_count1 = record1.observed_alleles
    allele_count2 = record2.observed_alleles
    if allele_count1 <= 1 or allele_count2 <= 1:
        return None

    alt_freq1 = record1.observed_alt / allele_count1
    alt_freq2 = record2.observed_alt / allele_count2
    numerator = ((alt_freq1 - alt_freq2) ** 2) - (
        (alt_freq1 * (1.0 - alt_freq1)) / (allele_count1 - 1)
    ) - ((alt_freq2 * (1.0 - alt_freq2)) / (allele_count2 - 1))
    denominator = (alt_freq1 * (1.0 - alt_freq2)) + (alt_freq2 * (1.0 - alt_freq1))
    return numerator, denominator


def calc_pi_windows(
    records: Iterable[Sequence[str] | str],
    window_size: int = 100000,
    populations: Optional[Sequence[str]] = None,
) -> list[dict[str, object]]:
    selected = set(populations) if populations else None
    windows: dict[tuple[str, str, int, int], dict[str, object]] = {}

    for locus in iter_loci(records):
        for record in locus:
            if selected is not None and record.population not in selected:
                continue
            if record.observed_alleles == 0:
                continue

            window_start, window_end = window_bounds(record.position, window_size)
            key = (record.population, record.chromosome, window_start, window_end)
            if key not in windows:
                windows[key] = {
                    "pop": record.population,
                    "chromosome": record.chromosome,
                    "window_pos_1": window_start,
                    "window_pos_2": window_end,
                    "pi": "NA",
                    "no_sites": 0,
                    "count_diffs": 0,
                    "count_comparisons": 0,
                    "count_missing": 0,
                }

            diffs, comparisons, missing, _, _ = _site_pi(record)
            window = windows[key]
            window["no_sites"] = int(window["no_sites"]) + 1
            window["count_diffs"] = int(window["count_diffs"]) + diffs
            window["count_comparisons"] = int(window["count_comparisons"]) + comparisons
            window["count_missing"] = int(window["count_missing"]) + missing

    output_rows = []
    for key in sorted(windows):
        row = windows[key]
        comparisons = int(row["count_comparisons"])
        if comparisons > 0:
            row["pi"] = int(row["count_diffs"]) / comparisons
        output_rows.append(row)

    return output_rows


def calc_tajima_d_windows(
    records: Iterable[Sequence[str] | str],
    window_size: int = 100000,
    populations: Optional[Sequence[str]] = None,
) -> list[dict[str, object]]:
    selected = set(populations) if populations else None
    windows: dict[tuple[str, str, int, int], dict[str, object]] = {}

    for locus in iter_loci(records):
        for record in locus:
            if selected is not None and record.population not in selected:
                continue

            window_start, window_end = window_bounds(record.position, window_size)
            key = (record.population, record.chromosome, window_start, window_end)
            if key not in windows:
                windows[key] = {
                    "pop": record.population,
                    "chromosome": record.chromosome,
                    "window_pos_1": window_start,
                    "window_pos_2": window_end,
                    "tajima_d": "NA",
                    "no_sites": 0,
                    "raw_pi": 0.0,
                    "raw_watterson_theta": 0.0,
                    "tajima_d_stdev": 0.0,
                    "_variant_sites": 0,
                    "_observed_allele_total_all": 0,
                    "_all_sites": 0,
                }

            diffs, comparisons, _, raw_pi, _ = _site_pi(record)
            window = windows[key]
            window["_all_sites"] = int(window["_all_sites"]) + 1
            window["_observed_allele_total_all"] = int(window["_observed_allele_total_all"]) + record.observed_alleles
            if record.observed_alleles > 0:
                window["no_sites"] = int(window["no_sites"]) + 1
            window["raw_pi"] = float(window["raw_pi"]) + raw_pi

            if record.observed_alt > 0:
                window["_variant_sites"] = int(window["_variant_sites"]) + 1
                window["raw_watterson_theta"] = float(window["raw_watterson_theta"]) + (1.0 / harmonic_number(record.observed_alleles))

    output_rows = []
    for key in sorted(windows):
        row = windows[key]
        all_sites = int(row["_all_sites"])
        no_sites = int(row["no_sites"])
        variant_sites = int(row["_variant_sites"])
        if all_sites > 0:
            mean_alleles = int(round(int(row["_observed_allele_total_all"]) / float(all_sites)))
        else:
            mean_alleles = 0

        if no_sites == 0 or mean_alleles < 2:
            row["tajima_d_stdev"] = "NA"

        if no_sites > 0 and mean_alleles >= 2 and variant_sites > 0:
            a1 = harmonic_number(mean_alleles)
            a2 = second_harmonic_number(mean_alleles)
            b1 = (mean_alleles + 1) / (3 * (mean_alleles - 1))
            b2 = (2 * (mean_alleles**2 + mean_alleles + 3)) / (9 * mean_alleles * (mean_alleles - 1))
            c1 = b1 - (1 / a1)
            c2 = b2 - ((mean_alleles + 2) / (a1 * mean_alleles)) + (a2 / (a1**2))
            e1 = c1 / a1
            e2 = c2 / (a1**2 + a2)
            d_stdev = sqrt((e1 * variant_sites) + (e2 * variant_sites * (variant_sites - 1)))
            row["tajima_d_stdev"] = d_stdev
            if d_stdev > 0:
                row["tajima_d"] = (
                    float(row["raw_pi"]) - float(row["raw_watterson_theta"])
                ) / d_stdev

        del row["_variant_sites"]
        del row["_observed_allele_total_all"]
        del row["_all_sites"]
        output_rows.append(row)

    return output_rows


def calc_dxy_windows(
    records: Iterable[Sequence[str] | str],
    window_size: int = 100000,
    populations: Optional[Sequence[str]] = None,
) -> list[dict[str, object]]:
    selected = set(populations) if populations else None
    windows: dict[tuple[str, str, str, int, int], dict[str, object]] = {}

    for locus in iter_loci(records):
        population_map = {record.population: record for record in locus}
        available_pops = sorted(population_map)
        if selected is not None:
            available_pops = [pop for pop in available_pops if pop in selected]

        for pop1, pop2 in combinations(available_pops, 2):
            record1 = population_map[pop1]
            record2 = population_map[pop2]
            diffs, comparisons, missing = _site_dxy(record1, record2)
            window_start, window_end = window_bounds(record1.position, window_size)
            key = (pop1, pop2, record1.chromosome, window_start, window_end)
            if key not in windows:
                windows[key] = {
                    "pop1": pop1,
                    "pop2": pop2,
                    "chromosome": record1.chromosome,
                    "window_pos_1": window_start,
                    "window_pos_2": window_end,
                    "avg_dxy": "NA",
                    "no_sites": 0,
                    "count_diffs": 0,
                    "count_comparisons": 0,
                    "count_missing": 0,
                }

            window = windows[key]
            if comparisons > 0:
                window["no_sites"] = int(window["no_sites"]) + 1
            window["count_diffs"] = int(window["count_diffs"]) + diffs
            window["count_comparisons"] = int(window["count_comparisons"]) + comparisons
            window["count_missing"] = int(window["count_missing"]) + missing

    output_rows = []
    for key in sorted(windows):
        row = windows[key]
        comparisons = int(row["count_comparisons"])
        if comparisons > 0:
            row["avg_dxy"] = int(row["count_diffs"]) / comparisons
        output_rows.append(row)

    return output_rows


def calc_dxy_windows_scantools_mean(
    records: Iterable[Sequence[str] | str],
    window_size: int = 100000,
    populations: Optional[Sequence[str]] = None,
) -> list[dict[str, object]]:
    selected = set(populations) if populations else None
    windows: dict[tuple[str, str, str, int, int], dict[str, object]] = {}

    for locus in iter_loci(records):
        population_map = {record.population: record for record in locus}
        available_pops = sorted(population_map)
        if selected is not None:
            available_pops = [pop for pop in available_pops if pop in selected]

        for pop1, pop2 in combinations(available_pops, 2):
            record1 = population_map[pop1]
            record2 = population_map[pop2]
            diffs, comparisons, missing = _site_dxy(record1, record2)
            window_start, window_end = window_bounds(record1.position, window_size)
            key = (pop1, pop2, record1.chromosome, window_start, window_end)
            if key not in windows:
                windows[key] = {
                    "pop1": pop1,
                    "pop2": pop2,
                    "chromosome": record1.chromosome,
                    "window_pos_1": window_start,
                    "window_pos_2": window_end,
                    "dxy_scantools_mean": "NA",
                    "no_sites": 0,
                    "sum_site_dxy": 0.0,
                    "count_diffs": 0,
                    "count_comparisons": 0,
                    "count_missing": 0,
                }

            window = windows[key]
            window["count_diffs"] = int(window["count_diffs"]) + diffs
            window["count_comparisons"] = int(window["count_comparisons"]) + comparisons
            window["count_missing"] = int(window["count_missing"]) + missing
            if comparisons > 0:
                window["no_sites"] = int(window["no_sites"]) + 1
                window["sum_site_dxy"] = float(window["sum_site_dxy"]) + (diffs / comparisons)

    output_rows = []
    for key in sorted(windows):
        row = windows[key]
        no_sites = int(row["no_sites"])
        if no_sites > 0:
            row["dxy_scantools_mean"] = float(row["sum_site_dxy"]) / no_sites
        output_rows.append(row)

    return output_rows


def calc_hudson_fst_windows(
    records: Iterable[Sequence[str] | str],
    window_size: int = 100000,
    populations: Optional[Sequence[str]] = None,
    minimum_snps: int = 1,
) -> list[dict[str, object]]:
    selected = set(populations) if populations else None
    windows: dict[tuple[str, str, str, int, int], dict[str, object]] = {}

    for locus in iter_loci(records):
        population_map = {record.population: record for record in locus}
        available_pops = sorted(population_map)
        if selected is not None:
            available_pops = [pop for pop in available_pops if pop in selected]

        for pop1, pop2 in combinations(available_pops, 2):
            site_components = _site_hudson_fst(population_map[pop1], population_map[pop2])
            if site_components is None:
                continue

            window_start, window_end = window_bounds(population_map[pop1].position, window_size)
            key = (pop1, pop2, population_map[pop1].chromosome, window_start, window_end)
            if key not in windows:
                windows[key] = {
                    "pop1": pop1,
                    "pop2": pop2,
                    "chromosome": population_map[pop1].chromosome,
                    "window_pos_1": window_start,
                    "window_pos_2": window_end,
                    "avg_hudson_fst": "NA",
                    "no_sites": 0,
                    "no_snps": 0,
                    "hudson_num": 0.0,
                    "hudson_den": 0.0,
                }

            numerator, denominator = site_components
            window = windows[key]
            window["no_sites"] = int(window["no_sites"]) + 1
            if denominator > 0.0:
                window["no_snps"] = int(window["no_snps"]) + 1
                window["hudson_num"] = float(window["hudson_num"]) + numerator
                window["hudson_den"] = float(window["hudson_den"]) + denominator

    output_rows = []
    for key in sorted(windows):
        row = windows[key]
        if int(row["no_snps"]) >= minimum_snps and float(row["hudson_den"]) > 0.0:
            row["avg_hudson_fst"] = float(row["hudson_num"]) / float(row["hudson_den"])
        output_rows.append(row)

    return output_rows