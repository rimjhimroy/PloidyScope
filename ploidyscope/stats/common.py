from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from typing import Iterable, Iterator, Optional, Sequence


MISSING_GENOTYPE = "-9"


@dataclass(frozen=True)
class PopulationRecord:
    population: str
    ploidy: int
    chromosome: str
    position: int
    genotypes: tuple[Optional[int], ...]
    observed_alt: int
    observed_alleles: int
    possible_alleles: int


def is_header(fields: Sequence[str]) -> bool:
    if not fields:
        return True
    first = fields[0].strip().lower()
    return first in {"pop", "population", "sample", "scaffold"}


def harmonic_number(n: int) -> float:
    if n <= 1:
        return 0.0
    return sum(1.0 / value for value in range(1, n))


def second_harmonic_number(n: int) -> float:
    if n <= 1:
        return 0.0
    return sum(1.0 / (value * value) for value in range(1, n))


def parse_record(fields: Sequence[str]) -> PopulationRecord:
    if len(fields) < 7:
        raise ValueError("Expected at least 7 tab-delimited columns per record")

    population = fields[0]
    ploidy = int(float(fields[1]))
    chromosome = fields[2]
    position = int(float(fields[3]))
    genotype_calls = []

    for value in fields[6:]:
        value = value.strip()
        if value == MISSING_GENOTYPE or value == "":
            genotype_calls.append(None)
        else:
            genotype_calls.append(int(float(value)))

    observed_calls = [value for value in genotype_calls if value is not None]
    observed_alt = sum(observed_calls)
    observed_alleles = len(observed_calls) * ploidy
    possible_alleles = len(genotype_calls) * ploidy

    return PopulationRecord(
        population=population,
        ploidy=ploidy,
        chromosome=chromosome,
        position=position,
        genotypes=tuple(genotype_calls),
        observed_alt=observed_alt,
        observed_alleles=observed_alleles,
        possible_alleles=possible_alleles,
    )


def build_record(
    population: str,
    ploidy: int,
    chromosome: str,
    position: int,
    genotype_calls: Sequence[Optional[int]],
) -> PopulationRecord:
    observed_calls = [value for value in genotype_calls if value is not None]
    observed_alt = sum(observed_calls)
    observed_alleles = len(observed_calls) * ploidy
    possible_alleles = len(genotype_calls) * ploidy

    return PopulationRecord(
        population=population,
        ploidy=ploidy,
        chromosome=chromosome,
        position=position,
        genotypes=tuple(genotype_calls),
        observed_alt=observed_alt,
        observed_alleles=observed_alleles,
        possible_alleles=possible_alleles,
    )


def load_population_map(path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(path, encoding="utf-8") as handle:
        sample_text = handle.read(2048)
        handle.seek(0)
        delimiter = "\t" if "\t" in sample_text else ","
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            if not row:
                continue
            lowered = [value.strip().lower() for value in row]
            if lowered[0] in {"sample", "sample_id", "individual"}:
                continue
            if len(row) < 2:
                raise ValueError("Population map must contain at least two columns: sample and population")
            sample = row[0].strip()
            population = row[1].strip()
            if not sample or not population:
                raise ValueError("Population map rows must have non-empty sample and population values")
            mapping[sample] = population
    if not mapping:
        raise ValueError("Population map is empty")
    return mapping


def iter_records(raw_records: Iterable[Sequence[str] | str]) -> Iterator[PopulationRecord]:
    for raw_record in raw_records:
        if isinstance(raw_record, PopulationRecord):
            yield raw_record
            continue
        if isinstance(raw_record, str):
            fields = raw_record.rstrip("\n").split("\t")
        else:
            fields = list(raw_record)

        if not fields or is_header(fields):
            continue

        yield parse_record(fields)


def iter_loci(raw_records: Iterable[Sequence[str] | str]) -> Iterator[list[PopulationRecord]]:
    current_key: Optional[tuple[str, int]] = None
    locus: list[PopulationRecord] = []

    for raw_record in raw_records:
        if isinstance(raw_record, list) and raw_record and isinstance(raw_record[0], PopulationRecord):
            if locus:
                yield locus
                locus = []
            yield raw_record
            current_key = None
            continue

        if isinstance(raw_record, tuple) and raw_record and isinstance(raw_record[0], PopulationRecord):
            if locus:
                yield locus
                locus = []
            yield list(raw_record)
            current_key = None
            continue

        for record in iter_records([raw_record]):
            key = (record.chromosome, record.position)
            if current_key is None:
                current_key = key

            if key != current_key:
                if locus:
                    yield locus
                locus = [record]
                current_key = key
                continue

            locus.append(record)

    if locus:
        yield locus


def iter_loci_from_vcf(
    vcf_path: str,
    population_map_path: str,
    populations: Optional[Sequence[str]] = None,
) -> Iterator[list[PopulationRecord]]:
    population_map = load_population_map(population_map_path)
    selected = set(populations) if populations else None
    sample_names: list[str]

    try:
        from cyvcf2 import VCF

        backend = "cyvcf2"
        vcf = VCF(vcf_path)
        sample_names = list(vcf.samples)
        variant_iterator = iter(vcf)
    except ModuleNotFoundError:
        try:
            from pysam import VariantFile
        except ModuleNotFoundError as exc:
            backend = "text"
            open_func = gzip.open if str(vcf_path).endswith(".gz") else open
            sample_names = []
            with open_func(vcf_path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        header_fields = line.rstrip("\n").split("\t")
                        sample_names = header_fields[9:]
                        break

            def _iter_text_variants():
                with open_func(vcf_path, "rt", encoding="utf-8") as handle:
                    for line in handle:
                        if line.startswith("#") or not line.strip():
                            continue
                        yield line.rstrip("\n").split("\t")

            variant_iterator = _iter_text_variants()
        else:
            backend = "pysam"
            vcf = VariantFile(vcf_path)
            sample_names = list(vcf.header.samples)
            variant_iterator = iter(vcf)

    sample_indices_by_population: dict[str, list[int]] = {}
    for index, sample in enumerate(sample_names):
        population = population_map.get(sample)
        if population is None:
            continue
        if selected is not None and population not in selected:
            continue
        sample_indices_by_population.setdefault(population, []).append(index)

    if not sample_indices_by_population:
        raise ValueError("No VCF samples matched the provided population map")

    for variant in variant_iterator:
        if backend == "cyvcf2":
            chromosome = variant.CHROM
            position = int(variant.POS)
            ref_allele = variant.REF
            alternates = variant.ALT
        elif backend == "pysam":
            chromosome = variant.contig
            position = int(variant.pos)
            ref_allele = variant.ref
            alternates = variant.alts
        else:
            chromosome = variant[0]
            position = int(variant[1])
            ref_allele = variant[3]
            alternates = () if variant[4] == "." else tuple(variant[4].split(","))
        if alternates is None:
            alternates = ()
        if len(alternates) > 1:
            continue
        if len(str(ref_allele)) != 1:
            continue
        if alternates and any(len(str(allele)) != 1 for allele in alternates):
            continue

        locus: list[PopulationRecord] = []
        for population in sorted(sample_indices_by_population):
            sample_indices = sample_indices_by_population[population]
            genotype_calls: list[Optional[int]] = []
            ploidy_values = set()

            for sample_index in sample_indices:
                if backend == "cyvcf2":
                    genotype = variant.genotypes[sample_index]
                    alleles = genotype[:-1]
                elif backend == "pysam":
                    sample_name = sample_names[sample_index]
                    genotype = variant.samples[sample_name].get("GT")
                    alleles = list(genotype) if genotype is not None else []
                else:
                    format_fields = variant[8].split(":")
                    sample_fields = variant[9 + sample_index].split(":")
                    format_map = dict(zip(format_fields, sample_fields))
                    gt_field = format_map.get("GT", ".")
                    separators_normalized = gt_field.replace("|", "/")
                    allele_fields = separators_normalized.split("/") if separators_normalized else []
                    alleles = []
                    for allele in allele_fields:
                        if allele == ".":
                            alleles.append(-1)
                        else:
                            alleles.append(int(allele))
                ploidy_values.add(len(alleles))
                if not alleles or any(allele is None or allele < 0 for allele in alleles):
                    genotype_calls.append(None)
                else:
                    genotype_calls.append(sum(1 for allele in alleles if allele == 1))

            nonzero_ploidies = {value for value in ploidy_values if value > 0}
            if not nonzero_ploidies:
                continue
            if len(nonzero_ploidies) > 1:
                raise ValueError(
                    f"Population {population} has inconsistent ploidy at {chromosome}:{position}"
                )

            ploidy = nonzero_ploidies.pop()
            locus.append(
                build_record(
                    population=population,
                    ploidy=ploidy,
                    chromosome=chromosome,
                    position=position,
                    genotype_calls=genotype_calls,
                )
            )

        if locus:
            yield locus


def window_bounds(position: int, window_size: int) -> tuple[int, int]:
    window_index = (position - 1) // window_size
    window_start = (window_index * window_size) + 1
    window_end = (window_index + 1) * window_size
    return window_start, window_end