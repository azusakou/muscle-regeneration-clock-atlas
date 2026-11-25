# Muscle Regeneration Clock Atlas

Expression-derived lineage phases, reference-relative alignment and programme-strength analyses for skeletal muscle repair.

The repository contains the analytical code, gene sets, study registry and compact inputs needed to reproduce the numerical sensitivity analyses and biological controls. Large expression matrices, manuscript files and reviewer correspondence are excluded. Raw data are available through the accessions in `config/datasets.tsv`.

## Requirements

R 4.4.2 and the packages recorded in `config/software_versions.tsv`. Core analyses use data.table, Matrix, SingleCellExperiment, scuttle, scDblFinder, scran, scater, mgcv and BiocParallel. Preprocessing also uses harmony, uwot, irlba and igraph; trajectory inference uses slingshot. Additional controls use edgeR, limma, rhdf5, DropletUtils and Seurat. Figures use ggplot2, patchwork and grid, with Quartz on macOS or Cairo elsewhere for embedded-font vector export. Package installation is not performed automatically.

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

A phase is an expression-state projection, not a longitudinal measurement of elapsed time. The ordered programmes overlap. Primary age comparisons use an injury-time anchor and exclude each evaluated young biological group from its own reference fit. The A/B and C/D geriatric capture pairs are conservative inferred groups, not verified animal identifiers; all three possible pairings are reported. Age-label tests refit the reference under each allocation. Within-accession lineage tests likewise refit all 15 allocations at days 1 and 3.5; an undefined allocation makes the corresponding test non-estimable. The depletion test also refits sham expectations and dispersion for all 20 three-of-six treatment allocations, retaining five eligible RSI animals. Pairing ranges are sensitivity ranges, not confidence intervals. Pooled libraries describe trajectories and do not supply animal replication.

The original benchmark uses unsupervised within-study programme standardization before accession-level splitting. It therefore evaluates normalized study transfer rather than fully inductive prediction. The central sensitivity analyses show that RSI depends materially on the phase estimator and reference. Depletion results are preliminary, and the negative whole-score contrast does not persist after removing the macrophage component. Gene coexpression does not demonstrate ligand secretion, receptor activity or causal cell communication.

## Independent biological controls

GSE180225 and GSE183802 are additional controls and do not enter the fixed atlas totals or original prediction benchmark. Download the GSE180225 RAW archive and unpack its 13 H5 files to `data/GSE180225/`. The author’s `Seurat_1.Rmd` in [Kurland_et_al_SCR_2023](https://github.com/jeku7901-CU/Kurland_et_al_SCR_2023) identifies `A_day4PI`/`A_day7PI` as young injured muscle and `OU_R1`/`OU_R2` as aged baseline, resolving conflicting original GEO titles. The mapping is explicit in the code.

```sh
Rscript R/biology.R . independent
```

The analysis reads sparse count arrays without loading millions of empty barcode strings, applies declared nucleus QC, removes doublets, scores the same programmes with Mymk/Tmem8c and Mymx/Gm7325 reconciliation, and fits simple linear young references for the three-time design. Young libraries are cross-fitted. A reference requires three observations at two times; every evaluated lineage requires 20 nuclei. Inadequate coverage suppresses RSI, including formal age tests with undefined permutation allocations. The separate myonuclear-phase test includes all 13 libraries and enumerates 360 within-time allocations.

For the denervation control, download the GSE183802 author Seurat RDS archive and set:

```sh
export MUSCLE_DENERVATION_RDS=/path/to/GSE183802_snRNA-seq_data_CTL_and_denervation.rds.gz
Rscript R/biology.R . denervation
```

The source has two gzip layers. The code retains author compartment labels, applies adaptive QC and sample-wise doublet removal, aggregates body myonuclei per animal, and runs edgeR/CAMERA. Fixed-four-lineage RSI is tested with all 20 reference-refitted condition allocations. An existing QC checkpoint can instead be selected through `MUSCLE_DENERVATION_QC`.

FAP response tests use the original six GSE291782 animals and the retained cell metadata generated by the `depletion` stage of `R/pipeline.R`:

```sh
export MUSCLE_DEPLETION_COUNTS=/path/to/GSE291782/extracted
export MUSCLE_DEPLETION_METADATA=results/GSE291782_failure/cell_metadata.tsv.gz
Rscript R/biology.R . depletion
```

All three FAP response panels belong to one multiple-testing family. A lineage-poor animal can contribute a FAP expression endpoint while remaining ineligible for RSI.

## Uniform-count NTX age sensitivity

Run the original preprocessing stages first to obtain the retained metadata for GSE159500 and GSE162172. Set the raw-count and metadata locations, then run:

```sh
export MUSCLE_DATA_ROOT=/path/to/muscle_regeneration_data
export MUSCLE_LIFESPAN_COUNTS=/path/to/GSE232106/extracted
export MUSCLE_NOTEXIN_COUNTS=/path/to/GSE143437
export MUSCLE_GSE159500_METADATA=results/GSE159500_extension/cell_metadata.tsv.gz
export MUSCLE_GSE162172_METADATA=results/GSE162172_old/cell_inclusion_and_scores.tsv.gz
Rscript R/biology.R . harmonized
Rscript R/biology.R . lineage_age
Rscript R/figures.R .
```

The same six programme-gene panels and log(1+CP10K) normalization are used across 65 NTX libraries. Official symbol synonyms are resolved before common coverage is required. Fifty-one conservative biological groups meet RSI criteria. A pooled-age reference cross-fits every evaluated group symmetrically and is independent of permuted age labels. Inference uses within-time permutations; lineage tests require two groups in every age category at an included time. `harmonized_scores` recomputes the final statistics from the compact programme summary after count extraction.

The denervation response is a positive biological control, while neither the RSI contrast in six denervation animals nor the age analyses establishes physiological validity. The retrieved human clinical supplement gives group means for the 72-subject bulk cohort, without individual functional values linked to the 17 nuclear donors; these summaries cannot be treated as donor-level outcomes. No new wet-lab experiments or matched individual force measurements were generated.
