# PloidyScope

PloidyScope is a small, deterministic Snakemake workflow and Python package for windowed population-genetic summaries on mixed-ploidy recoded genotype tables. It is designed as a straightforward implementation with one clean package path and one simple workflow entrypoint.

The codebase is intentionally direct:

- one package: `ploidyscope`
- one statistics namespace: `ploidyscope.stats`
- one table-driven Snakemake workflow
- explicit output tables with raw counts where relevant

## What is implemented

- `rho`: pairwise mixed-ploidy `rho` based directly on the original ScanTools `NestedAnova` logic from `pmonnahan/ScanTools`, with the same SPAGeDi-style variance decomposition.
- `dxy`: pairwise divergence with raw differences, raw comparisons, and missing-comparison counts, following the general output philosophy used by `pixy`.
- `pi`: within-population diversity with raw differences/comparisons/missing counts.
- `tajima_d`: within-population Tajima's D with explicit `raw_pi`, `raw_watterson_theta`, and `tajima_d_stdev` fields.

For Tajima's D, the denominator follows the post-issue-160 `pixy` direction: raw $\pi$ and raw Watterson's $\theta$ are kept in the numerator, while the denominator uses a single average observed allele count per window rather than the older per-allele-count summation.

## Input format

The workflow expects one tab-delimited table sorted by chromosome and position with this layout:

```text
pop    ploidy    chromosome    position    an    dp    genotype_1    genotype_2    ...
```

`genotype_*` values are per-individual alternate-allele counts. Missing genotypes must be encoded as `-9`, matching the original ScanTools recoded format.

This repository does not try to rebuild an end-to-end VCF preprocessing stack. The intended entrypoint is a recoded table produced upstream, either from existing preprocessing code or from your own table-generation step.

## Snakemake usage

Set `input_table` and `stats` in `config.yaml`, then run:

```bash
snakemake -n -p
snakemake --cores 4
```

By default the example config only requests `rho`, which keeps the workflow focused on the mixed-ploidy statistic most likely to need a dedicated run.

## CLI usage

You can run any single stat directly without Snakemake:

```bash
python -m ploidyscope.stats.cli \
  --stat rho \
  --infile data/example.recoded.tsv \
  --out results/ploidyscope/rho.tsv \
  --window-size 100000 \
  --minimum-snps 2
```

Supported `--stat` values are `rho`, `dxy`, `pi`, and `tajima_d`.

## Output tables

`rho.tsv`

```text
pop1    pop2    chromosome    window_pos_1    window_pos_2    rho    no_sites    no_snps    rho_num    rho_den    fst_num    fst_den
```

`dxy.tsv`

```text
pop1    pop2    chromosome    window_pos_1    window_pos_2    avg_dxy    no_sites    count_diffs    count_comparisons    count_missing
```

`pi.tsv`

```text
pop    chromosome    window_pos_1    window_pos_2    pi    no_sites    count_diffs    count_comparisons    count_missing
```

`tajima_d.tsv`

```text
pop    chromosome    window_pos_1    window_pos_2    tajima_d    no_sites    raw_pi    raw_watterson_theta    tajima_d_stdev
```

## Polyploid notes

- `rho` is the part of the codebase that stays closest to the original ScanTools variance decomposition, because that is the key mixed-ploidy statistic here.
- `dxy` and `pi` are computed from allele-count differences/comparisons, which extend naturally to arbitrary ploidy as long as the recoded table contains alternate-allele counts per individual.
- `tajima_d` is provided for polyploid recoded tables, but like the current `pixy` discussion, this remains a practical approximation under missing data rather than a final theoretical solution.

## Testing

From the repository root:

```bash
python -m pytest tests -q
```

## References

- Monnahan, P. et al. 2019. Pervasive population genomic consequences of genome duplication in Arabidopsis arenosa.
- ScanTools original implementation: `pmonnahan/ScanTools`
- `pixy` for output conventions and missing-data-aware raw summaries: `ksamuk/pixy`

