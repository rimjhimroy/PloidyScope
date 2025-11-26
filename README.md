# ScanTools Snakemake integration

This folder contains a Snakemake-ready wrapper around the original ScanTools scripts plus a refactored, testable subset of the code. The goal is to enable scalable, parallel execution of ScanTools-style windowed analyses while keeping a small, well-tested core that is easy to reuse in pipelines.

## What this repo contains
- `refactored/`: pure-function implementations of core algorithms (bpm, wpm) suitable for unit testing and Snakemake integration.
- `scripts/`: wrappers and the original ScanTools CLI for backward compatibility and quick CLI usage.
- `envs/scantools.yaml`: conda environment with pinned dependencies required to run the ScanTools python scripts.
- `config.yaml`: example configuration with paths to input VCF, reference, popkey file, and sliding-window settings.
- `tests/`: unit tests for the refactored modules (pytest).

## Design and workflow
1. Snakemake generates genomic windows (intervals) from the input VCF and launches per-window jobs using `bcftools` + the ScanTools wrapper or refactored CLI. This preserves parallelism and keeps the original ScanTools logic intact.
2. Each window is processed independently; results are written per-window to `results/` and later merged by Snakemake rules into final aggregates (windowed- and site-level outputs).
3. The refactored modules encapsulate the core algorithms so they can be unit-tested and imported in other scripts.

## Quick setup (conda)

1. Create and activate the environment:

```bash
conda env create -f envs/scantools.yaml
conda activate scantools
```

2. Or install with pip in an existing Python 3.9+ environment:

```bash
python -m pip install -r requirements.txt
```

## Example usage

Run the per-window pipeline (from repository root where Snakefile is present):

```bash
# dry-run to check DAG
snakemake -n -p

# run with 8 cores and create missing conda envs on-the-fly
snakemake -j 8 --use-conda
```

Run a single window with the refactored CLI (useful for debugging):

```bash
python3 scripts/run_scantools.py --vcf path/to.vcf.gz --chrom chr1 --start 1 --end 100000 --popkey configs/popkey.txt
```

## Testing

From this directory:

```bash
export PYTHONPATH=$(pwd)/../
python -m pytest tests -q
```

## Citation and References
This wrapper and the refactored code are intended as a convenience layer over the original ScanTools implementation. Please follow the license and citation requirements of the upstream project when redistributing or publishing results derived from this code.

If you use this Snakemake wrapper, please credit the wrapper author:

R.R. Choudhury (2025). ScanTools Snakemake wrapper. https://github.com/rimjhimroy/Scantools_snakemake

Suggested BibTeX for this wrapper:

```bibtex
@misc{rchoudhury_scantools_wrapper_2025,
	author = {Choudhury, R. R.},
	title = {ScanTools Snakemake wrapper},
	year = {2025},
	howpublished = {Repository / workflow in project},
	note = {URL: https://github.com/rimjhimroy/Scantools_snakemake}
}
```
Please credit the original ScanTools project when using this wrapper. The upstream project is hosted on GitHub:

Original ScanTools repository: https://github.com/pmonnahan/ScanTools (Paul Monnahan)

Monnahan, P. et al. Pervasive population genomic consequences of genome duplication in Arabidopsis arenosa. Nat. Ecol. Evol. 3, 457 (2019).

```bibtex
@article{monnahan_arenosa_2019,
	author = {Monnahan, P. and others},
	title = {Pervasive population genomic consequences of genome duplication in Arabidopsis arenosa},
	journal = {Nature Ecology & Evolution},
	year = {2019},
	volume = {3},
	pages = {457},
	note = {doi: 10.1038/s41559-019-0807-4}
}

```

## Notes and caveats
- The refactored modules are intended to be deterministic given fixed random seeds; set `numpy.random.seed()` in callers if exact reproducibility of downsampling is required.
- When running on real data, ensure sites excluded from analysis (low coverage, low mappability, repeats) are consistently masked when computing callable sites and summary statistics.

