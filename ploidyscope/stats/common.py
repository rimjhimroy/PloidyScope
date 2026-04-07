from __future__ import annotations

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


def iter_records(raw_records: Iterable[Sequence[str] | str]) -> Iterator[PopulationRecord]:
    for raw_record in raw_records:
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

    for record in iter_records(raw_records):
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


def window_bounds(position: int, window_size: int) -> tuple[int, int]:
    window_index = (position - 1) // window_size
    window_start = (window_index * window_size) + 1
    window_end = (window_index + 1) * window_size
    return window_start, window_end