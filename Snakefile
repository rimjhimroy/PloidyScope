configfile: "config.yaml"

from pathlib import Path


OUTDIR = Path(config.get("outdir", "results/ploidyscope"))
INPUT_TABLE = config.get("input_table")
INPUT_VCF = config.get("input_vcf")
POPULATION_MAP = config.get("population_map")
STATS = config.get("stats", ["rho"])
WINDOW_SIZE = int(config.get("window_size", 100000))
MINIMUM_SNPS = int(config.get("minimum_snps", 2))
POPULATIONS = config.get("populations", [])
CLI = "python -m ploidyscope.stats.cli"

if bool(INPUT_TABLE) == bool(INPUT_VCF):
    raise ValueError("Set exactly one of input_table or input_vcf in config.yaml")

if INPUT_VCF and not POPULATION_MAP:
    raise ValueError("population_map is required when input_vcf is set")

OUTDIR.mkdir(parents=True, exist_ok=True)

rule all:
    input:
        expand(str(OUTDIR / "{stat}.tsv"), stat=STATS)


rule compute_stat:
    input:
        lambda wildcards: [path for path in [INPUT_TABLE, INPUT_VCF, POPULATION_MAP] if path]
    output:
        str(OUTDIR / "{stat}.tsv")
    params:
        window_size=WINDOW_SIZE,
        minimum_snps=MINIMUM_SNPS,
        population_args=("--populations " + " ".join(POPULATIONS)) if POPULATIONS else "",
        input_args=(f"--infile {INPUT_TABLE}" if INPUT_TABLE else f"--vcf {INPUT_VCF} --popmap {POPULATION_MAP}"),
    wildcard_constraints:
        stat="|".join(STATS)
    shell:
        (
            "{CLI} --stat {wildcards.stat} {params.input_args} --out {output} "
            "--window-size {params.window_size} --minimum-snps {params.minimum_snps} "
            "{params.population_args}"
        )
