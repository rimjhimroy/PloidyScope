configfile: 'config.yaml'

import os
from pathlib import Path

OUTDIR = config['outdir']
VCF = config['vcf']
POPKEY = config['popkey']
REF = config['ref']
WINDOW = config.get('window_size', 100000)
STEP = config.get('step', WINDOW)
CHROMS = config.get('chroms', [])
SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'

os.makedirs(OUTDIR, exist_ok=True)
work_dir = Path('work')
work_dir.mkdir(exist_ok=True)

# helper: get intervals from VCF headers
def make_intervals_from_vcf(vcf, window, step):
    import subprocess
    import re

    # Try to get contig lengths from VCF header (preferred)
    try:
        header = subprocess.check_output(['bcftools', 'view', '-h', vcf], stderr=subprocess.DEVNULL).decode().splitlines()
    except Exception:
        header = []

    contigs = []
    contig_lengths = {}
    # parse lines like: ##contig=<ID=chr1,length=248956422>
    for line in header:
        line = line.strip()
        if not line.startswith('##contig=<'):
            continue
        inner = line[line.find('<')+1:line.rfind('>')]
        # split key=value pairs, allow commas inside quoted values conservatively
        parts = [p.strip() for p in inner.split(',') if p.strip()]
        name = None
        length = None
        for p in parts:
            if p.startswith('ID='):
                name = p.split('=',1)[1]
            elif p.startswith('length='):
                try:
                    length = int(p.split('=',1)[1])
                except Exception:
                    length = None
        if name:
            contigs.append(name)
            if length:
                contig_lengths[name] = length

    # If header had no contig lengths (or some contigs lack lengths), query VCF for observed max POS per contig
    try:
        out = subprocess.check_output(['bcftools', 'query', '-f', '%CHROM\t%POS\n', vcf], stderr=subprocess.DEVNULL).decode().splitlines()
        for ln in out:
            if not ln:
                continue
            chrom, pos = ln.split('\t')
            pos = int(pos)
            if chrom not in contigs:
                contigs.append(chrom)
            # record the maximum observed POS per contig
            contig_lengths[chrom] = max(contig_lengths.get(chrom, 0), pos)
    except Exception:
        # If bcftools query failed, try to at least get a unique list of chromosomes
        try:
            out = subprocess.check_output(['bcftools','query','-f','%CHROM\n', vcf], stderr=subprocess.DEVNULL).decode().splitlines()
            for chrom in out:
                if chrom not in contigs:
                    contigs.append(chrom)
        except Exception:
            # Nothing we can do; return an empty list
            return []

    intervals = []
    for chrom in contigs:
        length = contig_lengths.get(chrom)
        if not length:
            # If length unknown, create a single window starting at 1
            intervals.append((chrom, 1, window))
            continue
        # Create windows using 1-based coordinates: start..end inclusive
        start = 1
        while start <= length:
            end = min(start + window - 1, length)
            intervals.append((chrom, start, end))
            start += step

    return intervals

intervals = []
if len(CHROMS) > 0:
    # Prefer to generate full windows from the VCF and then filter to the
    # user-requested chromosomes so we get multiple windows per chromosome
    # where VCF contig lengths or positions are available.
    all_intervals = make_intervals_from_vcf(VCF, WINDOW, STEP)
    if all_intervals:
        # Filter generated intervals to only the requested chroms
        intervals = [iv for iv in all_intervals if iv[0] in set(CHROMS)]
    else:
        # Fallback to a single window per requested chromosome
        for chrom in CHROMS:
            intervals.append((chrom, 1, WINDOW))
else:
    intervals = make_intervals_from_vcf(VCF, WINDOW, STEP)

# Cache intervals to speed up snakemake parsing: write once and reuse
INTERVALS_CACHE = Path(OUTDIR) / 'intervals.tsv'
if INTERVALS_CACHE.exists():
    intervals = []
    for ln in INTERVALS_CACHE.read_text().splitlines():
        if not ln.strip():
            continue
        chrom, start, end = ln.split('\t')
        intervals.append((chrom, int(start), int(end)))
else:
    # write cache file so future snakemake parsing is fast
    try:
        INTERVALS_CACHE.parent.mkdir(parents=True, exist_ok=True)
        with INTERVALS_CACHE.open('w') as fh:
            for chrom, start, end in intervals:
                fh.write(f"{chrom}\t{start}\t{end}\n")
    except Exception:
        pass

# all window outputs for bpm
bpm_outs = [f"{OUTDIR}/bpm/{c[0]}_{c[1]}_{c[2]}_BPM.txt" for c in intervals]

USE_REFACTORED = config.get('use_refactored', False)

rule all:
    input:
        bpm_outs

rule window_vcf:
    output:
        vcf=temporary(OUTDIR + '/vcf_windows/{chr}_{start}_{end}.vcf.gz')
    params:
        outdir=OUTDIR
    conda:
        'envs/scantools.yaml'
    shell:
        'bcftools view -r {wildcards.chr}:{wildcards.start}-{wildcards.end} {VCF} | bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}'

rule run_scantools_per_window:
    input:
        vcf=lambda wildcards: OUTDIR + f"/vcf_windows/{wildcards.chr}_{wildcards.start}_{wildcards.end}.vcf.gz"
    output:
        OUTDIR + '/bpm/{chr}_{start}_{end}_BPM.txt'
    params:
        mode='bpm',
        outdir=OUTDIR
    conda:
        'envs/scantools.yaml'
    shell:
        'python {SCRIPTS}/run_scantools.py --mode {params.mode} --vcf {input.vcf} --popkey {POPKEY} --out {output} --window {wildcards.chr}:{wildcards.start}-{wildcards.end}'

if USE_REFACTORED:
    rule run_refactored_per_window:
        input:
            table=OUTDIR + '/tables/{chr}_{start}_{end}.table'
        output:
            OUTDIR + '/bpm/{chr}_{start}_{end}_BPM.txt'
        conda:
            'envs/scantools.yaml'
        shell:
            'python {SCRIPTS}/refactored/cli.py --mode bpm --infile {input.table} --out {output} --window {WINDOW}'
