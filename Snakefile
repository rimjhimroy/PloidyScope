configfile: "config.yaml"

from pathlib import Path


OUTDIR = Path(config.get("outdir", "results/ploidyscope"))
INPUT_TABLE = config["input_table"]
STATS = config.get("stats", ["rho"])
WINDOW_SIZE = int(config.get("window_size", 100000))
MINIMUM_SNPS = int(config.get("minimum_snps", 2))
POPULATIONS = config.get("populations", [])
CLI = "python -m ploidyscope.stats.cli"

OUTDIR.mkdir(parents=True, exist_ok=True)

rule all:
    input:
        expand(str(OUTDIR / "{stat}.tsv"), stat=STATS)


rule compute_stat:
    input:
        INPUT_TABLE
    output:
        str(OUTDIR / "{stat}.tsv")
    params:
        window_size=WINDOW_SIZE,
        minimum_snps=MINIMUM_SNPS,
        population_args=("--populations " + " ".join(POPULATIONS)) if POPULATIONS else "",
    wildcard_constraints:
        stat="|".join(STATS)
    shell:
        (
            "{CLI} --stat {wildcards.stat} --infile {input} --out {output} "
            "--window-size {params.window_size} --minimum-snps {params.minimum_snps} "
            "{params.population_args}"
        )
