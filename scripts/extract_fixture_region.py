#!/usr/bin/env python3
"""Extract a region from a VCF and write reusable fixture inputs.

Usage:
        python scripts/extract_fixture_region.py \
            --vcf data/input.vcf.gz \
            --metadata data/metadata.tsv \
            --region chr1:1-1000000 \
            --selected-pops POP1 POP2 \
            --outdir results/fixture_region
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import subprocess
from dataclasses import asdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_SOURCE_VCF = Path(
    "/data/users/rchoudhury/Biscutella_serpentine/results/gatk4/final_vcf/missingness_filtered.vcf.gz"
)
CANONICAL_SOURCE_METADATA = Path(
    "/data/users/rchoudhury/Biscutella_serpentine/configs/metadata_serpentine.tsv"
)
DEFAULT_COMPARISON_DIR = PROJECT_ROOT / "tests" / "data" / "fixture_comparison"
DEFAULT_FIXTURE_DIR = DEFAULT_COMPARISON_DIR / "fixture"
DEFAULT_EXTRACT_OUTDIR = PROJECT_ROOT / "results" / "fixture_region"
DEFAULT_OUTPUT_PREFIX = "fixture_region"
DEFAULT_REGION = "chr1:1-1000000"
DEFAULT_SITE_STRIDE = 50
CANONICAL_OUTPUT_PREFIX = "comparison_fixture"
CANONICAL_SELECTED_POPS = ("BL11-S2", "DKS", "TNA", "TDS")
CANONICAL_REGION = "Bv1:1-10000000"


@dataclass(frozen=True)
class SampleMetadata:
    sample: str
    population: str
    ploidy: int
    serpentine: str


@dataclass(frozen=True)
class PairPlan:
    pop1: str
    pop2: str
    ploidy1: int
    ploidy2: int
    route: str


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", default=str(DEFAULT_EXTRACT_OUTDIR))
    parser.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    parser.add_argument("--region", required=True)
    parser.add_argument("--site-stride", type=int, default=DEFAULT_SITE_STRIDE)
    parser.add_argument(
        "--index-python",
        default=str(PROJECT_ROOT / ".tool-baselines" / "pixy-venv" / "bin" / "python"),
    )
    parser.add_argument(
        "--selected-pops",
        nargs="+",
        required=True,
    )
    return parser


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=True)


def ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _get_required_value(row: dict[str, str], candidates: tuple[str, ...], label: str) -> str:
    for candidate in candidates:
        value = row.get(candidate)
        if value is not None and value.strip():
            return value.strip()
    raise ValueError(f"Metadata is missing a {label} column; tried {', '.join(candidates)}")


def parse_metadata(path: Path, selected_pops: set[str] | None = None) -> dict[str, SampleMetadata]:
    records: dict[str, SampleMetadata] = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Metadata file has no header: {path}")
        for raw_row in reader:
            row = {key.strip().lower(): value for key, value in raw_row.items() if key is not None}
            population = _get_required_value(row, ("pop", "population"), "population")
            if selected_pops is not None and population not in selected_pops:
                continue
            sample = _get_required_value(row, ("sample", "sample_id", "individual"), "sample")
            ploidy_text = _get_required_value(row, ("ploidy",), "ploidy")
            label = row.get("serpentine") or row.get("label") or row.get("group") or "unknown"
            record = SampleMetadata(
                sample=sample,
                population=population,
                ploidy=int(ploidy_text),
                serpentine=label.strip(),
            )
            existing = records.get(sample)
            if existing is not None and existing != record:
                raise ValueError(f"Inconsistent metadata for sample {sample}: {existing} vs {record}")
            records[sample] = record
    if not records:
        raise ValueError("No samples matched the requested populations")
    return records


def read_vcf_samples(vcf_path: Path) -> list[str]:
    with gzip.open(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"Could not find VCF header in {vcf_path}")


def region_chromosome(region: str) -> str:
    return region.split(":", 1)[0]


def parse_region(region: str) -> tuple[str, int, int]:
    chromosome, interval = region.split(":", 1)
    start_text, end_text = interval.split("-", 1)
    return chromosome, int(start_text), int(end_text)


def format_gt(call: tuple[int | None, ...] | None, phased: bool = False) -> str:
    if call is None:
        return "."
    separator = "|" if phased else "/"
    values = ["." if allele is None or allele < 0 else str(allele) for allele in call]
    return separator.join(values) if values else "."


def write_minimal_vcf_header(handle, chromosome: str, contig_length: int | None, sample_names: list[str]) -> None:
    handle.write("##fileformat=VCFv4.2\n")
    if contig_length is not None:
        handle.write(f"##contig=<ID={chromosome},length={contig_length}>\n")
    else:
        handle.write(f"##contig=<ID={chromosome}>\n")
    handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    if sample_names:
        handle.write("\t" + "\t".join(sample_names))
    handle.write("\n")


def write_minimal_vcf_record(handle, record, sample_names: list[str]) -> None:
    alt = "." if not record.alts else ",".join(record.alts)
    fields = [
        str(record.contig),
        str(record.pos),
        record.id if record.id is not None else ".",
        record.ref,
        alt,
        ".",
        "PASS",
        ".",
        "GT",
    ]
    sample_gts: list[str] = []
    for sample_name in sample_names:
        sample_data = record.samples[sample_name]
        gt = sample_data.get("GT")
        phased = bool(sample_data.phased) if hasattr(sample_data, "phased") else False
        sample_gts.append(format_gt(gt, phased=phased))
    handle.write("\t".join(fields + sample_gts) + "\n")


def open_vcf_writer(path: Path):
    if path.suffix != ".gz":
        return path.open("w", encoding="utf-8")
    from pysam import BGZFile

    raw_handle = BGZFile(str(path), "w")
    return io.TextIOWrapper(raw_handle, encoding="utf-8")


def extract_region_subset(
    input_vcf: Path,
    region: str,
    kept_samples: list[str],
    output_vcf: Path,
    site_stride: int,
) -> None:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    chromosome, start, end = parse_region(region)
    try:
        from pysam import VariantFile
    except ModuleNotFoundError:
        VariantFile = None

    if VariantFile is not None:
        with VariantFile(str(input_vcf)) as source:
            source.subset_samples(kept_samples)
            contig_length = None
            if chromosome in source.header.contigs:
                contig_length = source.header.contigs[chromosome].length
            with open_vcf_writer(output_vcf) as handle:
                write_minimal_vcf_header(handle, chromosome, contig_length, kept_samples)
                saw_variant = False
                for record_index, record in enumerate(source.fetch(chromosome, start - 1, end), start=1):
                    if record_index % site_stride != 0:
                        continue
                    write_minimal_vcf_record(handle, record, kept_samples)
                    saw_variant = True
        if not saw_variant:
            raise ValueError(f"No variants found in region {region}")
        return

    keep_indices: list[int] | None = None
    saw_variant = False
    in_target_chromosome = False
    with gzip.open(input_vcf, "rt", encoding="utf-8") as source, open_vcf_writer(output_vcf) as handle:
        record_index = 0
        for line in source:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                header_samples = fields[9:]
                keep_indices = [index for index, sample in enumerate(header_samples) if sample in kept_samples]
                write_minimal_vcf_header(handle, chromosome, None, [header_samples[index] for index in keep_indices])
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] != chromosome:
                if in_target_chromosome:
                    break
                continue
            in_target_chromosome = True
            position = int(fields[1])
            if position < start:
                continue
            if position > end:
                break
            record_index += 1
            if record_index % site_stride != 0:
                continue
            fixed = fields[:9]
            fixed[5] = "."
            fixed[6] = "PASS"
            fixed[7] = "."
            fixed[8] = "GT"
            sample_values = []
            for index in keep_indices:
                sample_field = fields[9 + index]
                sample_values.append(sample_field.split(":", 1)[0])
            handle.write("\t".join(fixed + sample_values) + "\n")
            saw_variant = True
    if keep_indices is None:
        raise ValueError(f"Could not read VCF header from {input_vcf}")
    if not saw_variant:
        raise ValueError(f"No variants found in region {region}")


def index_vcfgz(vcf_path: Path, python_exe: Path) -> None:
    code = f"import pysam; pysam.tabix_index(r'{vcf_path}', preset='vcf', force=True)"
    run([str(python_exe), "-c", code])


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line)
            handle.write("\n")


def write_metadata_subset(path: Path, samples: list[SampleMetadata]) -> None:
    rows = ["sample\tpop\tploidy\tserpentine"]
    rows.extend(
        f"{sample.sample}\t{sample.population}\t{sample.ploidy}\t{sample.serpentine}"
        for sample in samples
    )
    write_lines(path, rows)


def write_popmap(path: Path, samples: list[SampleMetadata]) -> None:
    rows = [f"{sample.sample}\t{sample.population}" for sample in samples]
    write_lines(path, rows)


def build_comparison_plan(samples: Iterable[SampleMetadata]) -> list[PairPlan]:
    ploidy_by_pop: dict[str, int] = {}
    for sample in samples:
        existing = ploidy_by_pop.get(sample.population)
        if existing is not None and existing != sample.ploidy:
            raise ValueError(f"Population {sample.population} has inconsistent ploidy")
        ploidy_by_pop[sample.population] = sample.ploidy

    plans: list[PairPlan] = []
    for pop1, pop2 in combinations(sorted(ploidy_by_pop), 2):
        ploidy1 = ploidy_by_pop[pop1]
        ploidy2 = ploidy_by_pop[pop2]
        route = "pixy" if ploidy1 == ploidy2 else "scantools"
        plans.append(PairPlan(pop1=pop1, pop2=pop2, ploidy1=ploidy1, ploidy2=ploidy2, route=route))
    return plans


def write_comparison_map(path: Path, plans: list[PairPlan]) -> None:
    rows = ["pop1\tpop2\tploidy1\tploidy2\troute"]
    rows.extend(
        f"{plan.pop1}\t{plan.pop2}\t{plan.ploidy1}\t{plan.ploidy2}\t{plan.route}"
        for plan in plans
    )
    write_lines(path, rows)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    input_vcf = Path(args.vcf)
    metadata_path = Path(args.metadata)
    outdir = Path(args.outdir)
    index_python = Path(args.index_python)
    ensure_exists(input_vcf, "input VCF")
    ensure_exists(metadata_path, "metadata")
    ensure_exists(index_python, "Python with pysam")

    selected_pops = set(args.selected_pops)
    metadata_records = parse_metadata(metadata_path, selected_pops)
    vcf_samples = set(read_vcf_samples(input_vcf))
    selected_samples = sorted(sample for sample in metadata_records if sample in vcf_samples)
    missing_samples = sorted(sample for sample in metadata_records if sample not in vcf_samples)
    if missing_samples:
        raise ValueError(f"Selected metadata samples missing from the VCF: {', '.join(missing_samples[:10])}")

    selected_metadata = [metadata_records[sample] for sample in selected_samples]
    comparison_plan = build_comparison_plan(selected_metadata)
    outdir.mkdir(parents=True, exist_ok=True)

    subset_vcfgz = outdir / f"{args.output_prefix}.vcf.gz"
    subset_metadata = outdir / f"{args.output_prefix}_metadata.tsv"
    subset_popmap = outdir / f"{args.output_prefix}_popmap.tsv"
    comparison_map = outdir / "comparison_map.tsv"

    extract_region_subset(input_vcf, args.region, selected_samples, subset_vcfgz, args.site_stride)
    index_vcfgz(subset_vcfgz, index_python)
    write_metadata_subset(subset_metadata, selected_metadata)
    write_popmap(subset_popmap, selected_metadata)
    write_comparison_map(comparison_map, comparison_plan)

    summary = {
        "fixture_vcf": str(subset_vcfgz),
        "metadata": str(subset_metadata),
        "popmap": str(subset_popmap),
        "comparison_map": str(comparison_map),
        "region": args.region,
        "site_stride": args.site_stride,
        "selected_populations": sorted(selected_pops),
        "plans": [asdict(plan) for plan in comparison_plan],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()