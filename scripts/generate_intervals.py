#!/usr/bin/env python3
"""Generate interval cache for ScanTools_snakemake Snakefile.

Usage: scripts/generate_intervals.py config.yaml
"""
import sys
from pathlib import Path
import yaml
import subprocess
import re


def make_intervals_from_vcf(vcf, window, step):
    try:
        header = (
            subprocess.check_output(
                ["bcftools", "view", "-h", vcf], stderr=subprocess.DEVNULL
            )
            .decode()
            .splitlines()
        )
    except Exception:
        header = []

    contigs = []
    contig_lengths = {}
    for line in header:
        line = line.strip()
        if not line.startswith("##contig=<"):
            continue
        inner = line[line.find("<") + 1 : line.rfind(">")]
        parts = [p.strip() for p in inner.split(",") if p.strip()]
        name = None
        length = None
        for p in parts:
            if p.startswith("ID="):
                name = p.split("=", 1)[1]
            elif p.startswith("length="):
                try:
                    length = int(p.split("=", 1)[1])
                except Exception:
                    length = None
        if name:
            contigs.append(name)
            if length:
                contig_lengths[name] = length

    try:
        out = (
            subprocess.check_output(
                ["bcftools", "query", "-f", "%CHROM\t%POS\n", vcf],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .splitlines()
        )
        for ln in out:
            if not ln:
                continue
            chrom, pos = ln.split("\t")
            pos = int(pos)
            if chrom not in contigs:
                contigs.append(chrom)
            contig_lengths[chrom] = max(contig_lengths.get(chrom, 0), pos)
    except Exception:
        try:
            out = (
                subprocess.check_output(
                    ["bcftools", "query", "-f", "%CHROM\n", vcf],
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .splitlines()
            )
            for chrom in out:
                if chrom not in contigs:
                    contigs.append(chrom)
        except Exception:
            return []

    intervals = []
    for chrom in contigs:
        length = contig_lengths.get(chrom)
        if not length:
            intervals.append((chrom, 1, window))
            continue
        start = 1
        while start <= length:
            end = min(start + window - 1, length)
            intervals.append((chrom, start, end))
            start += step
    return intervals


def main():
    cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
    vcf = cfg["vcf"]
    window = cfg.get("window_size", 100000)
    step = cfg.get("step", window)
    outdir = Path(cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    cache = outdir / "intervals.tsv"
    ints = make_intervals_from_vcf(vcf, window, step)
    with cache.open("w") as fh:
        for chrom, start, end in ints:
            fh.write(f"{chrom}\t{start}\t{end}\n")
    print("Wrote", cache, "with", len(ints), "intervals")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: generate_intervals.py config.yaml")
        sys.exit(1)
    main()
