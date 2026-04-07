from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Optional, Sequence

from .common import PopulationRecord
from .common import iter_loci
from .common import window_bounds


@dataclass(frozen=True)
class RhoSiteComponents:
    rho_num: float
    rho_den: float
    fst_num: float
    fst_den: float
    polymorphic: bool


def _filtered_genotypes(record: PopulationRecord) -> list[int]:
    return [value for value in record.genotypes if value is not None]


def calc_rho_site(records: Sequence[PopulationRecord]) -> Optional[RhoSiteComponents]:
    if len(records) < 2:
        return None

    locus = []
    for record in records:
        filtered = _filtered_genotypes(record)
        if not filtered:
            return None
        locus.append((record.population, record.ploidy, filtered))

    r = float(len(locus))
    p_i = []
    n_i = []
    ploidy_list = []
    p_ij = []
    ac_ij = []
    df_w = 0.0
    fs_nij2 = 0.0
    fs_snij2 = 0.0
    rs_snij2 = 0.0
    fs_snij2_over_snij = 0.0
    tac = 0
    tan = 0

    fssg = 0.0
    fssi = 0.0
    fssw = 0.0
    rssg = 0.0
    rssi = 0.0

    for _, ploidy, genotypes in locus:
        sample_size = float(len(genotypes))
        allele_number = float(ploidy * len(genotypes))
        alt_count = float(sum(genotypes))
        population_freq = alt_count / allele_number if allele_number else 0.0
        p_i.append(population_freq)
        n_i.append(sample_size)
        ploidy_list.append(float(ploidy))

        individual_freqs = []
        individual_counts = []
        fs_nij2_temp = 0.0
        population_alleles = 0

        for genotype in genotypes:
            individual_freq = float(genotype) / float(ploidy)
            individual_freqs.append(individual_freq)
            individual_counts.append(float(genotype))

            tac += int(genotype)
            tan += ploidy
            fs_nij2 += ploidy**2
            fs_nij2_temp += ploidy**2
            population_alleles += ploidy
            df_w += float(ploidy) - 1.0

        fs_snij2_over_snij += (
            fs_nij2_temp / population_alleles if population_alleles else 0.0
        )
        fs_snij2 += population_alleles**2
        rs_snij2 += sample_size**2
        p_ij.append(individual_freqs)
        ac_ij.append(individual_counts)

    if tan == 0:
        return None

    p_bar = float(tac) / float(tan)
    df_g = r - 1.0
    df_i = sum(value - 1.0 for value in n_i)
    if df_g <= 0.0 or df_i <= 0.0 or df_w <= 0.0:
        return None

    fn0bis = (fs_snij2_over_snij - (fs_nij2 / tan)) / df_g if tan else 0.0
    fn0 = (tan - fs_snij2_over_snij) / df_i
    fnb0 = (tan - (fs_snij2 / tan)) / df_g if tan else 0.0
    rnb0 = (sum(n_i) - (rs_snij2 / sum(n_i))) / df_g if sum(n_i) else 0.0

    if fn0 == 0.0 or fnb0 == 0.0 or rnb0 == 0.0:
        return None

    for population_index, individuals in enumerate(ac_ij):
        ploidy = int(ploidy_list[population_index])
        for individual_index, alt_count in enumerate(individuals):
            ref_count = ploidy - int(alt_count)
            for _ in range(ref_count):
                fssg += (p_i[population_index] - p_bar) ** 2
                fssi += (p_ij[population_index][individual_index] - p_i[population_index]) ** 2
                fssw += (0.0 - p_ij[population_index][individual_index]) ** 2
            for _ in range(int(alt_count)):
                fssg += (p_i[population_index] - p_bar) ** 2
                fssi += (p_ij[population_index][individual_index] - p_i[population_index]) ** 2
                fssw += (1.0 - p_ij[population_index][individual_index]) ** 2
            rssi += (p_ij[population_index][individual_index] - p_i[population_index]) ** 2
            rssg += (p_i[population_index] - p_bar) ** 2

    fms_g = fssg / df_g
    fms_i = fssi / df_i
    fms_w = fssw / df_w
    rms_g = rssg / df_g
    rms_i = rssi / df_i

    fs2_w = fms_w
    fs2_i = (fms_i - fs2_w) / fn0
    fs2_g = (fms_g - fs2_w - (fn0bis * fs2_i)) / fnb0
    rs2_i = rms_i
    rs2_g = (rms_g - rs2_i) / rnb0

    polymorphic = not (all(freq == 0.0 for freq in p_i) or all(freq == 1.0 for freq in p_i))

    return RhoSiteComponents(
        rho_num=rs2_g,
        rho_den=rs2_i + rs2_g,
        fst_num=fs2_g,
        fst_den=fs2_w + fs2_g + fs2_i,
        polymorphic=polymorphic,
    )


def calc_rho_windows(
    records: Iterable[Sequence[str] | str],
    window_size: int = 100000,
    minimum_snps: int = 2,
    populations: Optional[Sequence[str]] = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    selected = set(populations) if populations else None
    windows: dict[tuple[str, str, str, int, int], dict[str, object]] = {}
    genome = defaultdict(float)
    genome_sites: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)

    for locus in iter_loci(records):
        population_map = {record.population: record for record in locus}
        available_pops = sorted(population_map)
        if selected is not None:
            available_pops = [pop for pop in available_pops if pop in selected]

        for pop1, pop2 in combinations(available_pops, 2):
            pair_records = [population_map[pop1], population_map[pop2]]
            components = calc_rho_site(pair_records)
            if components is None:
                continue

            chromosome = pair_records[0].chromosome
            window_start, window_end = window_bounds(pair_records[0].position, window_size)
            key = (pop1, pop2, chromosome, window_start, window_end)
            if key not in windows:
                windows[key] = {
                    "pop1": pop1,
                    "pop2": pop2,
                    "chromosome": chromosome,
                    "window_pos_1": window_start,
                    "window_pos_2": window_end,
                    "rho": "NA",
                    "no_sites": 0,
                    "no_snps": 0,
                    "rho_num": 0.0,
                    "rho_den": 0.0,
                    "fst_num": 0.0,
                    "fst_den": 0.0,
                }

            window = windows[key]
            window["no_sites"] = int(window["no_sites"]) + 1
            genome_sites[(pop1, pop2)]["no_sites"] = genome_sites[(pop1, pop2)].get("no_sites", 0) + 1

            if components.polymorphic:
                window["no_snps"] = int(window["no_snps"]) + 1
                window["rho_num"] = float(window["rho_num"]) + components.rho_num
                window["rho_den"] = float(window["rho_den"]) + components.rho_den
                window["fst_num"] = float(window["fst_num"]) + components.fst_num
                window["fst_den"] = float(window["fst_den"]) + components.fst_den
                genome[(pop1, pop2, "rho_num")] += components.rho_num
                genome[(pop1, pop2, "rho_den")] += components.rho_den
                genome[(pop1, pop2, "fst_num")] += components.fst_num
                genome[(pop1, pop2, "fst_den")] += components.fst_den
                genome_sites[(pop1, pop2)]["no_snps"] = genome_sites[(pop1, pop2)].get("no_snps", 0) + 1

    output_rows = []
    for key in sorted(windows):
        row = windows[key]
        if int(row["no_snps"]) >= minimum_snps and float(row["rho_den"]) != 0.0:
            fac = float(row["rho_num"]) / float(row["rho_den"])
            row["rho"] = fac / (1.0 + fac)
        output_rows.append(row)

    genomewide = []
    for pop1, pop2 in sorted({(key[0], key[1]) for key in genome_sites}):
        rho_den = genome[(pop1, pop2, "rho_den")]
        fst_den = genome[(pop1, pop2, "fst_den")]
        fac = genome[(pop1, pop2, "rho_num")] / rho_den if rho_den else None
        genomewide.append(
            {
                "pop1": pop1,
                "pop2": pop2,
                "rho": (fac / (1.0 + fac)) if fac is not None else "NA",
                "rho_num": genome[(pop1, pop2, "rho_num")],
                "rho_den": rho_den,
                "fst": (
                    genome[(pop1, pop2, "fst_num")] / fst_den if fst_den else "NA"
                ),
                "fst_num": genome[(pop1, pop2, "fst_num")],
                "fst_den": fst_den,
                "no_sites": genome_sites[(pop1, pop2)].get("no_sites", 0),
                "no_snps": genome_sites[(pop1, pop2)].get("no_snps", 0),
            }
        )

    return output_rows, {"pairs": genomewide}