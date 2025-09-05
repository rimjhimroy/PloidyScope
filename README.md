ScanTools Snakemake integration

This folder contains a Snakemake-ready wrapper around the original ScanTools scripts plus a refactored, testable subset of the code.

Layout
- `refactored/`: pure-function implementations of core algorithms (bpm, wpm) suitable for unit testing and Snakemake integration.
- `scripts/`: wrappers and the original ScanTools CLI for backward compatibility.
- `envs/scantools.yaml`: conda environment with dependencies needed to run the ScanTools python scripts (numpy pinned).
- `config.yaml`: example config with paths to input VCF, reference, popkey file, and window settings.
- `tests/`: unit tests for the refactored modules.

How it works
1. Snakemake splits the input VCF into genomic windows using bcftools view with intervals defined from a sliding-window generator in the Snakefile. This preserves parallelism and avoids modifying the original ScanTools splitting logic.
2. Each window can be processed independently using either the original wrapper (`scripts/run_scantools.py`) or the new refactored CLI.
3. Outputs are produced per-window and later concatenated/merged by Snakemake rules into final results.

Quick setup (conda)

1. Create the environment:

```bash
conda env create -f envs/scantools.yaml
conda activate scantools
```

2. Alternatively, use pip in an existing Python 3.9+ environment:

```bash
python -m pip install -r requirements.txt
```

Running tests

From this directory run:

```bash
export PYTHONPATH=$(pwd)/../
python -m pytest tests -q
```

Notes

- The tests add a small sys.path alias to handle case-sensitive package directory names (`ScanTools_snakemake` vs `scantools_snakemake`).
- The refactored modules are intended to be deterministic given fixed random seeds; set `numpy.random.seed()` in callers if exact reproducibility of downsampling is required.
