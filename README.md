# Muscle Regeneration Clock Atlas

Expression-derived lineage phases, reference-relative alignment and programme-strength analyses for skeletal muscle repair.

The repository contains the analytical code, gene sets, study registry and compact inputs needed to reproduce the numerical sensitivity analyses and biological controls. Large expression matrices, manuscript files and reviewer correspondence are excluded. Raw data are available through the accessions in `config/datasets.tsv`.

## Requirements

R 4.4.2 and the packages recorded in `config/software_versions.tsv`. Core analyses use data.table, Matrix, SingleCellExperiment, scuttle, scDblFinder, scran, scater, mgcv and BiocParallel. Preprocessing also uses harmony, uwot, irlba and igraph; trajectory inference uses slingshot. Figures use ggplot2, patchwork and grid, with Quartz on macOS or Cairo elsewhere for embedded-font vector export. Package installation is not performed automatically.

## Reproduce the sensitivity analyses

Run from the repository root:

```sh
Rscript R/sensitivity.R .
```

The first command refits independent injury, consensus and leave-one-lineage-out references; compares parameters and programme omissions; enumerates age and depletion permutations; and performs 100 equal-depth resamples. Results are written to `results/`. Numerical table inputs retain full precision. The compressed cell-metadata object contains only retained annotations and the six programme scores required for cell-recovery resampling; it contains no expression-count matrix.

## Reproduce the gene-level controls

Download and unpack GSE232106 author-deposited SoupX-corrected matrices. Set the directory containing the geriatric matrix bundles and `young_soupx` subdirectory:

```sh
export MUSCLE_LIFESPAN_COUNTS=/path/to/GSE232106/extracted
Rscript R/biology.R .
Rscript R/figures.R .
```

The analysis matches the deposited cell identifiers to retained singlets, uses log(1 + 10,000 * count / library total), and evaluates inflammatory, restorative, FAP, Schwann and postsynaptic NMJ modules. Per-sample gene tables are checkpoints. The manuscript figure directory can be selected with `MUSCLE_FIGURE_DIR`; otherwise figures are written inside the repository.

## Rebuild from raw data

`R/pipeline.R` contains the original ingestion, annotation and modelling stages in one file. List stages with `Rscript R/pipeline.R help`. Set `MUSCLE_DATA_ROOT` for downloaded/extracted data and `MUSCLE_PROCESSED_ROOT` for processed objects. Core stages are:

1. `metadata INPUT.soft.gz OUTPUT.tsv` to parse sample metadata.
2. `prepare_cardiotoxin`, `prepare_notexin`, `prepare_lifespan` and `prepare_volumetric` to ingest the core studies; explicit input and output paths can be supplied after each stage name.
3. `annotate INPUT_QC.rds OUTPUT_ANNOTATED.rds` for each core study.
4. `integrate config/core_state_space_inputs.tsv OUTPUT.rds` to construct the common state space.
5. `discovery`, `multiview`, `depletion`, `hypothyroid`, `notexin_extension`, `barium_chloride`, `old_notexin`, `human_exercise`, `spatial`, `human_age`, `dynamics`, `programmes`, `curio` and `trajectory` to run the retained analyses.

The metadata stage writes the sample tables expected in `results/00_registry/`, including the GSE143437 and GSE215808 sample tables. Independent cohort bundles belong in `MUSCLE_DATA_ROOT/extracted/ACCESSION`; the human-age stage uses `data/GSE167186/` as in its original implementation. The spatial stage uses verified coordinate files; expression-only GSE266933 objects are excluded from neighbourhood analysis. Source matrices are deliberately not mirrored here.

## Interpretation

A phase is an expression-state projection, not a longitudinal measurement of elapsed time. The ordered programmes overlap. Primary age comparisons use an injury-time anchor and exclude each evaluated young biological group from its own reference fit. The A/B and C/D geriatric capture pairs are conservative inferred groups, not verified animal identifiers. Pooled libraries describe trajectories and do not supply animal replication.

The original benchmark uses unsupervised within-study programme standardization before accession-level splitting. It therefore evaluates normalized study transfer rather than fully inductive prediction. The central sensitivity analyses show that RSI depends materially on the phase estimator and reference. Depletion results are preliminary, and the negative whole-score contrast does not persist after removing the macrophage component. Gene coexpression does not demonstrate ligand secretion, receptor activity or causal cell communication.
