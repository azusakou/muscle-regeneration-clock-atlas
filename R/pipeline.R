#!/usr/bin/env Rscript
# Muscle Regeneration Clock Atlas: preprocessing and analytical stages.
argv <- base::commandArgs(trailingOnly=TRUE)
stage <- if(length(argv)) argv[1] else "help"
stage_args <- if(length(argv)>1) argv[-1] else character()
commandArgs <- function(trailingOnly=FALSE) if(trailingOnly) stage_args else base::commandArgs()
options(stringsAsFactors=FALSE)
# Sample-level clock metrics. Cells are inputs to phase estimation, never replicates.

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2L) return(NA_real_)
  mu <- weighted_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])
}

# deviation = observed cell-type phase minus the young-reference phase expected
# at the sample's consensus biological time. This preserves normal lineage lags.
regenerative_synchrony_index <- function(deviation, reliability = NULL,
                                         sigma_reference = 0.15,
                                         min_cell_types = 4L) {
  if (is.null(reliability)) reliability <- rep(1, length(deviation))
  ok <- is.finite(deviation) & is.finite(reliability) & reliability > 0
  if (sum(ok) < min_cell_types) return(NA_real_)
  dispersion <- weighted_var(deviation[ok], reliability[ok])
  exp(-dispersion / sigma_reference^2)
}

phase_abnormalities <- function(deviation, persistence = NULL,
                                delay_cut = 0.15, advance_cut = -0.15,
                                persistence_cut = 0.20) {
  cls <- rep("on_time", length(deviation))
  cls[is.finite(deviation) & deviation >= delay_cut] <- "phase_delay"
  cls[is.finite(deviation) & deviation <= advance_cut] <- "phase_advance"
  if (!is.null(persistence)) {
    cls[is.finite(persistence) & persistence >= persistence_cut] <- "temporal_persistence"
  }
  cls[!is.finite(deviation)] <- NA_character_
  cls
}

# Pairwise lag is directional and defined on an open regeneration interval.
# Do not wrap it modulo 2*pi: early inflammation and late remodeling are not adjacent.
pairwise_phase_lag <- function(phase_a, phase_b) phase_a - phase_b

assert_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Missing required R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

parse_injury_time_days <- function(sample_id) {
  key <- trimws(as.character(sample_id))
  out <- rep(NA_real_, length(key))
  out[grepl("non.?inj|uninj|control|d0$", key, ignore.case = TRUE)] <- 0

  numeric_token <- sub(
    ".*?([0-9]+(?:[.-][0-9]+)?).*",
    "\\1",
    key
  )
  numeric_token <- gsub("-", ".", numeric_token, fixed = TRUE)
  numeric_time <- suppressWarnings(as.numeric(numeric_token))
  out[is.na(out)] <- numeric_time[is.na(out)]
  out[grepl("hour|hpi|hours", key, ignore.case = TRUE)] <-
    out[grepl("hour|hpi|hours", key, ignore.case = TRUE)] / 24
  out
}

robust_qc_flags <- function(sce, sample_id) {
  stopifnot(length(sample_id) == ncol(sce))
  batch <- factor(sample_id)

  low_sum <- scuttle::isOutlier(
    sce$sum,
    batch = batch,
    type = "lower",
    log = TRUE,
    nmads = 3
  )
  low_detected <- scuttle::isOutlier(
    sce$detected,
    batch = batch,
    type = "lower",
    log = TRUE,
    nmads = 3
  )
  high_sum <- scuttle::isOutlier(
    sce$sum,
    batch = batch,
    type = "higher",
    log = TRUE,
    nmads = 5
  )
  high_detected <- scuttle::isOutlier(
    sce$detected,
    batch = batch,
    type = "higher",
    log = TRUE,
    nmads = 5
  )
  high_mt <- scuttle::isOutlier(
    sce$subsets_Mt_percent,
    batch = batch,
    type = "higher",
    log = FALSE,
    nmads = 3
  )

  sce$qc_low_sum <- low_sum
  sce$qc_low_detected <- low_detected
  sce$qc_high_sum <- high_sum
  sce$qc_high_detected <- high_detected
  sce$qc_high_mt <- high_mt
  sce$passes_core_qc <- !(
    low_sum | low_detected | high_sum | high_detected | high_mt
  )

  thresholds <- list(
    low_sum = attr(low_sum, "thresholds"),
    low_detected = attr(low_detected, "thresholds"),
    high_sum = attr(high_sum, "thresholds"),
    high_detected = attr(high_detected, "thresholds"),
    high_mt = attr(high_mt, "thresholds")
  )
  list(sce = sce, thresholds = thresholds)
}

sparse_from_dataframe <- function(x, chunk_size = 1000L) {
  stopifnot(is.data.frame(x), chunk_size > 0)
  starts <- seq.int(1L, ncol(x), by = chunk_size)
  chunks <- vector("list", length(starts))

  for (i in seq_along(starts)) {
    first <- starts[[i]]
    last <- min(first + chunk_size - 1L, ncol(x))
    message(
      "Converting columns ", first, "-", last, " of ", ncol(x),
      " to sparse"
    )
    dense_chunk <- as.matrix(x[, first:last, drop = FALSE])
    storage.mode(dense_chunk) <- "double"
    chunks[[i]] <- methods::as(dense_chunk, "dgCMatrix")
    rm(dense_chunk)
    invisible(gc(FALSE))
  }

  out <- do.call(cbind, chunks)
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  out
}

write_session_manifest <- function(path, extra = list()) {
  manifest <- c(
    list(
      generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      r_version = R.version.string,
      platform = R.version$platform,
      session_info = utils::capture.output(utils::sessionInfo())
    ),
    extra
  )
  jsonlite::write_json(
    manifest,
    path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
}

read_marker_table <- function(path, name_col, marker_col) {
  x <- data.table::fread(path, sep = "\t", data.table = FALSE)
  stopifnot(name_col %in% names(x), marker_col %in% names(x))
  sets <- strsplit(x[[marker_col]], ";", fixed = TRUE)
  names(sets) <- x[[name_col]]
  lapply(sets, function(z) unique(trimws(z[nzchar(trimws(z))])))
}

robust_z <- function(x) {
  center <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 1e-8) {
    spread <- stats::sd(x, na.rm = TRUE)
  }
  if (!is.finite(spread) || spread <= 1e-8) return(rep(0, length(x)))
  (x - center) / spread
}

bounded_standard_z <- function(x, limit = 5) {
  center <- mean(x, na.rm = TRUE)
  spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread <= 1e-8) return(rep(0, length(x)))
  out <- (x - center) / spread
  pmax(-limit, pmin(limit, out))
}

score_gene_sets <- function(log_expr, gene_sets) {
  genes <- rownames(log_expr)
  out <- matrix(
    NA_real_,
    nrow = ncol(log_expr),
    ncol = length(gene_sets),
    dimnames = list(colnames(log_expr), names(gene_sets))
  )
  coverage <- data.frame(
    program = names(gene_sets),
    requested = lengths(gene_sets),
    present = 0L,
    present_genes = NA_character_
  )

  for (i in seq_along(gene_sets)) {
    present <- intersect(gene_sets[[i]], genes)
    coverage$present[[i]] <- length(present)
    coverage$present_genes[[i]] <- paste(present, collapse = ";")
    if (length(present)) {
      out[, i] <- Matrix::colMeans(log_expr[present, , drop = FALSE])
    }
  }
  list(scores = out, coverage = coverage)
}

assign_cluster_labels <- function(score_z, clusters, min_top = 0, min_margin = 0.15) {
  dt <- data.table::as.data.table(score_z)
  dt[, cluster_id := as.character(clusters)]
  cluster_scores <- dt[, lapply(.SD, stats::median, na.rm = TRUE),
                       by = cluster_id,
                       .SDcols = colnames(score_z)]

  score_matrix <- as.matrix(cluster_scores[, -1, with = FALSE])
  top_index <- max.col(score_matrix, ties.method = "first")
  top_score <- score_matrix[cbind(seq_len(nrow(score_matrix)), top_index)]
  second_score <- apply(score_matrix, 1, function(z) sort(z, decreasing = TRUE)[2])
  margin <- top_score - second_score
  label <- colnames(score_matrix)[top_index]
  label[top_score < min_top | margin < min_margin] <- "Unknown_or_ambiguous"

  cluster_map <- data.frame(
    cluster_id = cluster_scores$cluster_id,
    coarse_cell_type = label,
    annotation_top_score = top_score,
    annotation_margin = margin,
    annotation_confidence = stats::plogis(margin),
    stringsAsFactors = FALSE
  )
  cell_map <- cluster_map[match(as.character(clusters), cluster_map$cluster_id), ]
  list(cluster_map = cluster_map, cell_map = cell_map, cluster_scores = cluster_scores)
}


stages <- list(
  "metadata" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript scripts/04_parse_geo_soft.R input.soft.gz output.tsv",
    call. = FALSE
  )
}
input <- args[[1]]
output <- args[[2]]
if (!file.exists(input)) stop("Missing input: ", input)

lines <- readLines(gzfile(input), warn = FALSE)
starts <- grep("^\\^SAMPLE = ", lines)
if (!length(starts)) stop("No SAMPLE records in ", input)
ends <- c(starts[-1L] - 1L, length(lines))

clean_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

records <- lapply(seq_along(starts), function(i) {
  block <- lines[starts[[i]]:ends[[i]]]
  accession <- sub("^\\^SAMPLE = ", "", block[[1]])
  out <- list(sample_accession = accession)

  scalar_fields <- c(
    title = "^!Sample_title = ",
    source_name = "^!Sample_source_name_ch1 = ",
    organism = "^!Sample_organism_ch1 = ",
    library_strategy = "^!Sample_library_strategy = ",
    library_source = "^!Sample_library_source = ",
    library_selection = "^!Sample_library_selection = "
  )
  for (field in names(scalar_fields)) {
    hit <- grep(scalar_fields[[field]], block, value = TRUE)
    out[[field]] <- if (length(hit)) {
      sub(scalar_fields[[field]], "", hit[[1]])
    } else {
      NA_character_
    }
  }

  characteristics <- sub(
    "^!Sample_characteristics_ch1 = ",
    "",
    grep("^!Sample_characteristics_ch1 = ", block, value = TRUE)
  )
  for (value in characteristics) {
    pieces <- strsplit(value, ":", fixed = TRUE)[[1]]
    if (length(pieces) >= 2L) {
      key <- clean_key(pieces[[1]])
      val <- trimws(paste(pieces[-1L], collapse = ":"))
      if (!is.null(out[[key]]) && !is.na(out[[key]])) {
        out[[key]] <- paste(out[[key]], val, sep = " | ")
      } else {
        out[[key]] <- val
      }
    }
  }

  supplements <- sub(
    "^!Sample_supplementary_file_[0-9]+ = ",
    "",
    grep("^!Sample_supplementary_file_[0-9]+ = ", block, value = TRUE)
  )
  out$supplementary_files <- if (length(supplements)) {
    paste(supplements, collapse = " | ")
  } else {
    NA_character_
  }
  out
})

all_names <- unique(unlist(lapply(records, names)))
records <- lapply(records, function(x) {
  missing <- setdiff(all_names, names(x))
  x[missing] <- NA_character_
  as.data.frame(x[all_names], check.names = FALSE)
})
samples <- data.table::rbindlist(records, fill = TRUE)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(samples, output, sep = "\t", na = "NA")
cat("Parsed", nrow(samples), "GEO sample records to", output, "\n")

  },
  "prepare_cardiotoxin" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
data_root <- Sys.getenv("MUSCLE_DATA_ROOT", "data")
processed_root <- Sys.getenv(
  "MUSCLE_PROCESSED_ROOT", file.path(data_root, "processed")
)
input <- if (length(args) >= 1L) args[[1]] else
  file.path(data_root, "raw/GSE138826/GSE138826_regen_data.rds.gz")
out_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(processed_root, "GSE138826")
chunk_size <- if (length(args) >= 3L) as.integer(args[[3]]) else 750L

assert_packages(c(
  "R.utils", "Biobase", "Matrix", "SingleCellExperiment", "scuttle",
  "scDblFinder", "BiocParallel", "data.table", "jsonlite"
))

if (!file.exists(input)) stop("Missing input: ", input)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/01_ingest/GSE138826", recursive = TRUE, showWarnings = FALSE)

set.seed(18018)
tmp <- tempfile(fileext = ".rds")
on.exit(unlink(tmp), add = TRUE)
R.utils::gunzip(input, destname = tmp, overwrite = TRUE, remove = FALSE)

message("Reading deposited GSE138826 object")
deposited <- readRDS(tmp)
if (!is.list(deposited) || !"RNA" %in% names(deposited)) {
  stop("Unexpected GSE138826 object structure")
}

rna <- deposited[["RNA"]]
rm(deposited)
invisible(gc())

counts_df <- Biobase::pData(rna)
if (!is.data.frame(counts_df)) stop("RNA payload is not a data.frame")
if (anyNA(counts_df)) stop("RNA count payload contains NA values")

counts <- sparse_from_dataframe(counts_df, chunk_size = chunk_size)
rm(counts_df, rna)
invisible(gc())

sample_id <- sub("_.*$", "", colnames(counts))
sample_id <- trimws(sample_id)
time_days <- parse_injury_time_days(sample_id)

col_data <- S4Vectors::DataFrame(
  dataset_id = "GSE138826",
  cell_id = colnames(counts),
  sample_id = sample_id,
  capture_id = sample_id,
  injury_model = "CTX",
  time_label = sample_id,
  time_days = time_days,
  species = "Mus musculus",
  assay_type = "scRNA-seq",
  biological_unit_note = "Each time-point library pools three male mice"
)
rownames(col_data) <- colnames(counts)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts),
  colData = col_data
)
rm(counts)
invisible(gc())

mt <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
sce <- scuttle::addPerCellQC(sce, subsets = list(Mt = mt))

qc <- robust_qc_flags(sce, sce$sample_id)
sce <- qc$sce

message("Running scDblFinder by capture library")
sce <- scDblFinder::scDblFinder(
  sce,
  samples = "capture_id",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = TRUE
)
sce$passes_qc_and_singlet <- sce$passes_core_qc &
  sce$scDblFinder.class == "singlet"

saveRDS(sce, file.path(out_dir, "GSE138826_qc_flagged_sce.rds"), compress = FALSE)

summary_dt <- data.table::as.data.table(
  as.data.frame(SummarizedExperiment::colData(sce))
)[, .(
  input_cells = .N,
  low_sum = sum(qc_low_sum),
  low_detected = sum(qc_low_detected),
  high_sum = sum(qc_high_sum),
  high_detected = sum(qc_high_detected),
  high_mt = sum(qc_high_mt),
  core_qc_pass = sum(passes_core_qc),
  doublets = sum(scDblFinder.class == "doublet"),
  final_singlet_qc_pass = sum(passes_qc_and_singlet),
  median_umis = as.numeric(stats::median(sum)),
  median_genes = as.numeric(stats::median(detected)),
  median_mt_percent = as.numeric(stats::median(subsets_Mt_percent))
), by = .(dataset_id, sample_id, time_days)]

data.table::fwrite(
  summary_dt,
  "results/01_ingest/GSE138826/qc_summary_by_sample.tsv",
  sep = "\t",
  na = "NA"
)
saveRDS(
  qc$thresholds,
  "results/01_ingest/GSE138826/qc_thresholds.rds",
  compress = "xz"
)

metadata_dt <- data.table::as.data.table(
  as.data.frame(SummarizedExperiment::colData(sce)),
  keep.rownames = "cell_barcode"
)
data.table::fwrite(
  metadata_dt,
  "results/01_ingest/GSE138826/cell_metadata.tsv.gz",
  sep = "\t",
  na = "NA"
)

write_session_manifest(
  "results/01_ingest/GSE138826/analysis_manifest.json",
  list(
    dataset_id = "GSE138826",
    input = normalizePath(input),
    output = normalizePath(
      file.path(out_dir, "GSE138826_qc_flagged_sce.rds"),
      mustWork = FALSE
    ),
    chunk_size = chunk_size,
    qc_policy = paste(
      "Per-library robust outliers: 3 MAD lower for log total counts and",
      "detected genes; 5 MAD upper for log total counts and detected genes;",
      "3 MAD upper for mitochondrial percentage."
    ),
    doublet_method = "scDblFinder per capture library",
    ambient_rna = paste(
      "Not used as a hard filter because raw empty-droplet matrices were not",
      "included in the deposited processed object."
    )
  )
)

message("Completed GSE138826 ingest and QC flagging")

  },
  "prepare_notexin" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
data_root <- Sys.getenv("MUSCLE_DATA_ROOT", "data")
processed_root <- Sys.getenv(
  "MUSCLE_PROCESSED_ROOT", file.path(data_root, "processed")
)
input_dir <- if (length(args) >= 1L) args[[1]] else
  file.path(data_root, "raw/GSE143437")
out_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(processed_root, "GSE143437")
chunk_size <- if (length(args) >= 3L) as.integer(args[[3]]) else 500L

assert_packages(c(
  "Matrix", "SingleCellExperiment", "scuttle", "scDblFinder",
  "BiocParallel", "data.table", "jsonlite"
))

count_file <- file.path(
  input_dir,
  "GSE143437_DeMicheli_MuSCatlas_rawdata.txt.gz"
)
cell_metadata_file <- file.path(
  input_dir,
  "GSE143437_DeMicheli_MuSCatlas_metadata.txt.gz"
)
sample_metadata_file <- "results/00_registry/GSE143437_samples.tsv"
needed <- c(count_file, cell_metadata_file, sample_metadata_file)
if (any(!file.exists(needed))) {
  stop("Missing inputs: ", paste(needed[!file.exists(needed)], collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- "results/01_ingest/GSE143437"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(18018)

message("Reading GSE143437 wide raw-count table")
count_dt <- data.table::fread(
  count_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)
gene <- make.unique(as.character(count_dt[[1]]))
count_dt[[1]] <- NULL
count_df <- as.data.frame(count_dt, check.names = FALSE)
rm(count_dt)
rownames(count_df) <- gene
counts <- sparse_from_dataframe(count_df, chunk_size = chunk_size)
rm(count_df)
invisible(gc())

cell_meta <- data.table::fread(
  cell_metadata_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  data.table = TRUE
)
data.table::setnames(cell_meta, 1L, "cell_id")
if (!setequal(colnames(counts), cell_meta$cell_id)) {
  stop("Count-matrix cells do not match deposited metadata")
}
cell_meta <- cell_meta[match(colnames(counts), cell_id)]

sample_meta <- data.table::fread(
  sample_metadata_file,
  sep = "\t",
  data.table = TRUE
)
sample_meta <- sample_meta[, .(
  sample_id = title,
  sample_accession,
  sex,
  age,
  genotype,
  geo_injury = injury
)]
cell_meta <- merge(
  cell_meta,
  sample_meta,
  by.x = "sampleID",
  by.y = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
cell_meta <- cell_meta[match(colnames(counts), cell_id)]
if (anyNA(cell_meta$sample_accession)) {
  stop("Failed to map one or more cells to GEO sample metadata")
}

time_days <- suppressWarnings(
  as.numeric(sub(".*?([0-9]+).*", "\\1", cell_meta$geo_injury))
)
time_days[
  grepl("No injury|Day 0", cell_meta$geo_injury, ignore.case = TRUE)
] <- 0

col_data <- S4Vectors::DataFrame(
  dataset_id = "GSE143437",
  cell_id = cell_meta$cell_id,
  sample_id = cell_meta$sampleID,
  sample_accession = cell_meta$sample_accession,
  capture_id = cell_meta$sampleID,
  injury_model = "notexin",
  time_label = cell_meta$geo_injury,
  time_days = time_days,
  sex = cell_meta$sex,
  age = cell_meta$age,
  genotype = cell_meta$genotype,
  author_cell_annotation = cell_meta$cell_annotation,
  author_nUMI = cell_meta$nUMI,
  author_nGene = cell_meta$nGene,
  author_percent_mito_fraction = cell_meta$percent_mito,
  species = "Mus musculus",
  assay_type = "scRNA-seq",
  biological_unit_note = "One mouse-derived capture library per GEO sample"
)
rownames(col_data) <- colnames(counts)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts),
  colData = col_data
)
rm(counts)
invisible(gc())

mt <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
sce <- scuttle::addPerCellQC(sce, subsets = list(Mt = mt))
qc <- robust_qc_flags(sce, sce$sample_id)
sce <- qc$sce

sce <- scDblFinder::scDblFinder(
  sce,
  samples = "capture_id",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = TRUE
)
sce$passes_qc_and_singlet <- sce$passes_core_qc &
  sce$scDblFinder.class == "singlet"

saveRDS(
  sce,
  file.path(out_dir, "GSE143437_qc_flagged_sce.rds"),
  compress = FALSE
)

metadata <- data.table::as.data.table(
  as.data.frame(SummarizedExperiment::colData(sce)),
  keep.rownames = "cell_barcode"
)
summary <- metadata[, .(
  input_cells = .N,
  low_sum = sum(qc_low_sum),
  low_detected = sum(qc_low_detected),
  high_sum = sum(qc_high_sum),
  high_detected = sum(qc_high_detected),
  high_mt = sum(qc_high_mt),
  core_qc_pass = sum(passes_core_qc),
  doublets = sum(scDblFinder.class == "doublet"),
  final_singlet_qc_pass = sum(passes_qc_and_singlet),
  median_umis = as.numeric(stats::median(sum)),
  median_genes = as.numeric(stats::median(detected)),
  median_mt_percent = as.numeric(stats::median(subsets_Mt_percent))
), by = .(
  dataset_id, sample_id, sample_accession, time_days, sex, age
)]

data.table::fwrite(
  summary,
  file.path(result_dir, "qc_summary_by_sample.tsv"),
  sep = "\t"
)
data.table::fwrite(
  metadata,
  file.path(result_dir, "cell_metadata.tsv.gz"),
  sep = "\t"
)
saveRDS(
  qc$thresholds,
  file.path(result_dir, "qc_thresholds.rds"),
  compress = "xz"
)

write_session_manifest(
  file.path(result_dir, "analysis_manifest.json"),
  list(
    dataset_id = "GSE143437",
    inputs = normalizePath(needed),
    output = normalizePath(
      file.path(out_dir, "GSE143437_qc_flagged_sce.rds"),
      mustWork = FALSE
    ),
    chunk_size = chunk_size,
    qc_policy = paste(
      "Per-library robust outliers: 3 MAD lower for log total counts and",
      "detected genes; 5 MAD upper for log total counts and detected genes;",
      "3 MAD upper for mitochondrial percentage."
    ),
    doublet_method = "scDblFinder per capture library",
    annotation_source = "Deposited author annotation retained for validation",
    design_caveat = paste(
      "D0 n=3, D2 n=2, D5 n=3, D7 n=2 with partial sex/age imbalance;",
      "all models must include available covariates and sensitivity analyses."
    )
  )
)
message("Completed GSE143437 ingest and QC flagging")

  },
  "prepare_lifespan" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
legacy_root <- Sys.getenv("MUSCLE_LEGACY_DATA_ROOT", "data")
data_root <- Sys.getenv("MUSCLE_DATA_ROOT", "data")
processed_root <- Sys.getenv(
  "MUSCLE_PROCESSED_ROOT", file.path(data_root, "processed")
)
input_dir <- if (length(args) >= 1L) args[[1]] else
  Sys.getenv("MUSCLE_LIFESPAN_COUNTS", file.path(legacy_root, "raw/mouse_lifespan_regeneration/extracted"))
out_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(processed_root, "GSE232106")

assert_packages(c(
  "Matrix", "SingleCellExperiment", "scuttle", "scDblFinder",
  "BiocParallel", "data.table", "jsonlite"
))
if (!dir.exists(input_dir)) stop("Missing input directory: ", input_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- "results/01_ingest/GSE232106"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(18018)

read_mtx_bundle <- function(matrix_file, feature_file, barcode_file, sample_id) {
  message("Reading ", sample_id)
  con <- if (grepl("\\.gz$", matrix_file)) gzfile(matrix_file, "rt") else matrix_file
  x <- Matrix::readMM(con)
  x <- methods::as(x, "dgCMatrix")
  feature <- data.table::fread(
    feature_file, header = FALSE, data.table = FALSE, showProgress = FALSE
  )
  barcode <- data.table::fread(
    barcode_file, header = FALSE, data.table = FALSE, showProgress = FALSE
  )
  gene_col <- if (ncol(feature) >= 2L) 2L else 1L
  gene <- make.unique(as.character(feature[[gene_col]]))
  if (nrow(x) != length(gene) || ncol(x) != nrow(barcode)) {
    stop("Dimension mismatch in ", sample_id)
  }
  rownames(x) <- gene
  colnames(x) <- paste(sample_id, barcode[[1]], sep = "__")
  x
}

ger_matrix <- Sys.glob(file.path(input_dir, "GSM*_soupx_matrix.mtx.gz"))
ger_matrix <- ger_matrix[!grepl("/\\._", ger_matrix)]
ger_prefix <- sub("_soupx_matrix\\.mtx\\.gz$", "", ger_matrix)
ger_ids <- basename(ger_prefix)

young_dirs <- Sys.glob(file.path(input_dir, "young_soupx", "GSM*"))
young_dirs <- young_dirs[dir.exists(young_dirs) & !grepl("/\\._", young_dirs)]
young_ids <- basename(young_dirs)

if (length(ger_ids) != 24L || length(young_ids) != 8L) {
  stop(
    "Expected 24 geriatric and 8 young libraries; found ",
    length(ger_ids), " and ", length(young_ids)
  )
}

matrices <- vector("list", length(ger_ids) + length(young_ids))
names(matrices) <- c(ger_ids, young_ids)
for (i in seq_along(ger_ids)) {
  matrices[[ger_ids[[i]]]] <- read_mtx_bundle(
    paste0(ger_prefix[[i]], "_soupx_matrix.mtx.gz"),
    paste0(ger_prefix[[i]], "_soupx_genes.tsv.gz"),
    paste0(ger_prefix[[i]], "_soupx_barcodes.tsv.gz"),
    ger_ids[[i]]
  )
}
for (i in seq_along(young_ids)) {
  d <- young_dirs[[i]]
  matrices[[young_ids[[i]]]] <- read_mtx_bundle(
    file.path(d, "matrix.mtx"),
    file.path(d, "genes.tsv"),
    file.path(d, "barcodes.tsv"),
    young_ids[[i]]
  )
}

common_genes <- Reduce(intersect, lapply(matrices, rownames))
if (length(common_genes) < 10000L) {
  stop("Unexpectedly small shared gene universe: ", length(common_genes))
}
matrices <- lapply(matrices, function(x) x[common_genes, , drop = FALSE])
counts <- do.call(cbind, matrices)
rm(matrices)
invisible(gc())

sample_id <- sub("__.*$", "", colnames(counts))
age_group <- ifelse(grepl("_Ger_", sample_id), "geriatric", "young")
time_token <- sub(
  ".*_(D[0-9]+(?:_[0-9]+)?).*", "\\1", sample_id, perl = TRUE
)
time_days <- as.numeric(sub("_", ".", sub("^D", "", time_token)))
replicate <- sub(".*_", "", sample_id)

col_data <- S4Vectors::DataFrame(
  dataset_id = "GSE232106",
  cell_id = colnames(counts),
  sample_id = sample_id,
  sample_accession = sub("_.*$", "", sample_id),
  capture_id = sample_id,
  injury_model = "notexin",
  time_label = time_token,
  time_days = time_days,
  age_group = age_group,
  replicate = replicate,
  species = "Mus musculus",
  assay_type = "scRNA-seq",
  count_layer = "author-deposited SoupX-corrected counts",
  biological_unit_note = "One mouse-derived capture library per GEO sample"
)
rownames(col_data) <- colnames(counts)

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts),
  colData = col_data
)
rm(counts)
invisible(gc())

mt <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
sce <- scuttle::addPerCellQC(sce, subsets = list(Mt = mt))
qc <- robust_qc_flags(sce, sce$sample_id)
sce <- qc$sce

message("Running scDblFinder independently within 32 capture libraries")
sce <- scDblFinder::scDblFinder(
  sce,
  samples = "capture_id",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = TRUE
)
sce$passes_qc_and_singlet <- sce$passes_core_qc &
  sce$scDblFinder.class == "singlet"

output <- file.path(out_dir, "GSE232106_qc_flagged_sce.rds")
base::saveRDS(sce, output, compress = FALSE)

metadata <- data.table::as.data.table(
  as.data.frame(SummarizedExperiment::colData(sce)),
  keep.rownames = "cell_barcode"
)
summary <- metadata[, .(
  input_cells = .N,
  low_sum = sum(qc_low_sum),
  low_detected = sum(qc_low_detected),
  high_sum = sum(qc_high_sum),
  high_detected = sum(qc_high_detected),
  high_mt = sum(qc_high_mt),
  core_qc_pass = sum(passes_core_qc),
  doublets = sum(scDblFinder.class == "doublet"),
  final_singlet_qc_pass = sum(passes_qc_and_singlet),
  median_umis = as.numeric(stats::median(sum)),
  median_genes = as.numeric(stats::median(detected)),
  median_mt_percent = as.numeric(stats::median(subsets_Mt_percent))
), by = .(
  dataset_id, sample_id, sample_accession, age_group, time_days, replicate
)]

data.table::fwrite(
  summary, file.path(result_dir, "qc_summary_by_sample.tsv"), sep = "\t"
)
data.table::fwrite(
  metadata, file.path(result_dir, "cell_metadata.tsv.gz"), sep = "\t"
)
base::saveRDS(
  qc$thresholds, file.path(result_dir, "qc_thresholds.rds"), compress = "xz"
)

write_session_manifest(
  file.path(result_dir, "analysis_manifest.json"),
  list(
    dataset_id = "GSE232106",
    input = normalizePath(input_dir),
    output = normalizePath(output, mustWork = FALSE),
    libraries = length(unique(sample_id)),
    input_cells = ncol(sce),
    shared_genes = length(common_genes),
    count_layer = "author-deposited SoupX-corrected matrices",
    qc_policy = paste(
      "Per-library robust outliers: 3 MAD lower for log total counts and",
      "detected genes; 5 MAD upper for log total counts and detected genes;",
      "3 MAD upper for mitochondrial percentage."
    ),
    doublet_method = "scDblFinder per capture library",
    design_note = paste(
      "Geriatric: four mice at D0, D1, D2, D3.5, D5 and D7;",
      "young: four mice at D1 and D3.5."
    )
  )
)
message("Completed GSE232106 ingest and QC flagging")

  },
  "prepare_volumetric" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
data_root <- Sys.getenv("MUSCLE_DATA_ROOT", "data")
processed_root <- Sys.getenv(
  "MUSCLE_PROCESSED_ROOT", file.path(data_root, "processed")
)
input_dir <- if (length(args) >= 1L) args[[1]] else
  file.path(data_root, "raw/GSE215808")
out_dir <- if (length(args) >= 2L) args[[2]] else
  file.path(processed_root, "GSE215808")

assert_packages(c(
  "DropletUtils", "SingleCellExperiment", "scuttle", "scDblFinder",
  "BiocParallel", "data.table", "jsonlite"
))

files <- c(
  GSM8043019 = file.path(
    input_dir,
    "GSM8043019_7698-JL-1-filtered_feature_bc_matrix.h5"
  ),
  GSM8043020 = file.path(
    input_dir,
    "GSM8043020_7698-JL-2-filtered_feature_bc_matrix.h5"
  )
)
if (any(!file.exists(files))) {
  stop("Missing GSE215808 H5 inputs: ", paste(files[!file.exists(files)], collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
result_dir <- "results/01_ingest/GSE215808"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(18018)

condition_map <- c(
  GSM8043019 = "VML_vehicle",
  GSM8043020 = "VML_Maresin1"
)

sce <- DropletUtils::read10xCounts(
  samples = files,
  sample.names = names(files),
  col.names = TRUE,
  row.names = "symbol",
  type = "HDF5",
  version = "3"
)
sce$dataset_id <- "GSE215808"
sce$cell_id <- colnames(sce)
sce$sample_id <- as.character(sce$Sample)
sce$capture_id <- sce$sample_id
sce$condition <- unname(condition_map[sce$sample_id])
sce$injury_model <- "VML"
sce$time_label <- "Day7"
sce$time_days <- 7
sce$species <- "Mus musculus"
sce$assay_type <- "scRNA-seq"
sce$biological_unit_note <- paste(
  "One deposited library per condition; no replicate-level differential",
  "inference is permitted within this dataset"
)

mt <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
sce <- scuttle::addPerCellQC(sce, subsets = list(Mt = mt))
qc <- robust_qc_flags(sce, sce$sample_id)
sce <- qc$sce

sce <- scDblFinder::scDblFinder(
  sce,
  samples = "capture_id",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = TRUE
)
sce$passes_qc_and_singlet <- sce$passes_core_qc &
  sce$scDblFinder.class == "singlet"

saveRDS(
  sce,
  file.path(out_dir, "GSE215808_qc_flagged_sce.rds"),
  compress = FALSE
)

metadata <- data.table::as.data.table(
  as.data.frame(SummarizedExperiment::colData(sce)),
  keep.rownames = "cell_barcode"
)
summary <- metadata[, .(
  input_cells = .N,
  low_sum = sum(qc_low_sum),
  low_detected = sum(qc_low_detected),
  high_sum = sum(qc_high_sum),
  high_detected = sum(qc_high_detected),
  high_mt = sum(qc_high_mt),
  core_qc_pass = sum(passes_core_qc),
  doublets = sum(scDblFinder.class == "doublet"),
  final_singlet_qc_pass = sum(passes_qc_and_singlet),
  median_umis = as.numeric(stats::median(sum)),
  median_genes = as.numeric(stats::median(detected)),
  median_mt_percent = as.numeric(stats::median(subsets_Mt_percent))
), by = .(dataset_id, sample_id, condition, time_days)]

data.table::fwrite(
  summary,
  file.path(result_dir, "qc_summary_by_sample.tsv"),
  sep = "\t"
)
data.table::fwrite(
  metadata,
  file.path(result_dir, "cell_metadata.tsv.gz"),
  sep = "\t"
)
saveRDS(
  qc$thresholds,
  file.path(result_dir, "qc_thresholds.rds"),
  compress = "xz"
)
write_session_manifest(
  file.path(result_dir, "analysis_manifest.json"),
  list(
    dataset_id = "GSE215808",
    inputs = normalizePath(files),
    output = normalizePath(
      file.path(out_dir, "GSE215808_qc_flagged_sce.rds"),
      mustWork = FALSE
    ),
    qc_policy = paste(
      "Per-library robust outliers: 3 MAD lower for log total counts and",
      "detected genes; 5 MAD upper for log total counts and detected genes;",
      "3 MAD upper for mitochondrial percentage."
    ),
    doublet_method = "scDblFinder per deposited capture library",
    ambient_rna = paste(
      "Not used as a hard filter because raw empty-droplet matrices were not",
      "deposited."
    ),
    inference_constraint = paste(
      "The two conditions lack deposited independent biological replicates;",
      "use for state/perturbation concordance only, not stand-alone p-values."
    )
  )
)

message("Completed GSE215808 ingest and QC flagging")

  },
  "annotate" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste(
      "Usage: Rscript scripts/05_annotate_core_datasets.R",
      "input_qc_sce.rds output_annotated_sce.rds"
    ),
    call. = FALSE
  )
}
input <- args[[1]]
output <- args[[2]]
dataset_id <- sub("_qc_flagged_sce\\.rds$", "", basename(input))
result_dir <- file.path("results/02_annotation", dataset_id)

assert_packages(c(
  "SingleCellExperiment", "scuttle", "scran", "scater", "Matrix",
  "igraph", "data.table", "jsonlite", "BiocParallel", "BiocSingular"
))
if (!file.exists(input)) stop("Missing input: ", input)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(18018)

sce <- readRDS(input)
required_metadata <- c("passes_qc_and_singlet", "sample_id")
if (!all(required_metadata %in% colnames(SummarizedExperiment::colData(sce)))) {
  stop("Input lacks QC metadata: ", paste(required_metadata, collapse = ", "))
}
sce <- sce[, sce$passes_qc_and_singlet]
if (!ncol(sce)) stop("No QC-passing singlets")

message("Log-normalizing ", format(ncol(sce), big.mark = ","), " cells")
sce <- scuttle::logNormCounts(sce)

message("Selecting batch-aware highly variable genes")
var_fit <- scran::modelGeneVar(
  sce,
  block = factor(sce$sample_id),
  assay.type = "logcounts",
  BPPARAM = BiocParallel::SerialParam()
)
nuisance <- grepl(
  "^(mt-|Rpl|Rps|Hba-|Hbb-|Gm[0-9]+)",
  rownames(var_fit),
  ignore.case = TRUE
)
ordered <- order(var_fit$bio, decreasing = TRUE, na.last = NA)
hvg <- rownames(var_fit)[ordered]
hvg <- hvg[!nuisance[match(hvg, rownames(var_fit))]]
hvg <- head(hvg, 3000L)
if (length(hvg) < 500L) stop("Too few HVGs after nuisance exclusion")

message("Running PCA; sample is not regressed because it is confounded with time")
sce <- scater::runPCA(
  sce,
  subset_row = hvg,
  ncomponents = 50,
  BSPARAM = BiocSingular::IrlbaParam(),
  BPPARAM = BiocParallel::SerialParam()
)

graph <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 30)
cluster <- igraph::cluster_louvain(graph)$membership
sce$cluster_id <- factor(cluster)

sce <- scater::runUMAP(
  sce,
  dimred = "PCA",
  ncomponents = 2,
  n_neighbors = 30,
  min_dist = 0.3,
  BPPARAM = BiocParallel::SerialParam()
)

cell_sets <- read_marker_table(
  "config/cell_type_marker_sets.tsv",
  "cell_type",
  "positive_markers"
)
cell_negative <- read_marker_table(
  "config/cell_type_marker_sets.tsv",
  "cell_type",
  "negative_or_exclusion_markers"
)
cell_scored <- score_gene_sets(SummarizedExperiment::assay(sce, "logcounts"), cell_sets)
negative_scored <- score_gene_sets(
  SummarizedExperiment::assay(sce, "logcounts"),
  cell_negative
)
cell_score <- cell_scored$scores - 0.5 * negative_scored$scores
cell_score_z <- apply(cell_score, 2, robust_z)
rownames(cell_score_z) <- colnames(sce)

labels <- assign_cluster_labels(cell_score_z, sce$cluster_id)
sce$coarse_cell_type <- labels$cell_map$coarse_cell_type
sce$annotation_top_score <- labels$cell_map$annotation_top_score
sce$annotation_margin <- labels$cell_map$annotation_margin
sce$annotation_confidence <- labels$cell_map$annotation_confidence
sce$annotation_backend <- "conservative cluster-level canonical-marker scoring"

state_sets <- read_marker_table(
  "config/state_gene_sets.tsv",
  "state_name",
  "positive_markers"
)
state_scored <- score_gene_sets(
  SummarizedExperiment::assay(sce, "logcounts"),
  state_sets
)
state_z <- apply(state_scored$scores, 2, robust_z)
rownames(state_z) <- colnames(sce)
for (state in colnames(state_z)) {
  SummarizedExperiment::colData(sce)[[paste0("state_score_", state)]] <-
    state_z[, state]
}

# These are evidence axes only. Final state labels require time, lineage,
# persistence, spatial or perturbation evidence in later stages.
sce$provisional_state_axis <- colnames(state_z)[max.col(state_z)]
sce$state_assignment_status <- "provisional_expression_axis_not_fate"

# Counts remain in the immutable QC checkpoint. The annotated checkpoint stores
# logcounts plus metadata/reductions only, avoiding a second multi-gigabyte copy
# of the raw count assay. base::saveRDS also preserves delayed/HDF5 references.
sce_out <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    logcounts = SummarizedExperiment::assay(sce, "logcounts")
  ),
  rowData = SummarizedExperiment::rowData(sce),
  colData = SummarizedExperiment::colData(sce),
  reducedDims = SingleCellExperiment::reducedDims(sce)
)
base::saveRDS(sce_out, output, compress = FALSE)
rm(sce_out)

data.table::fwrite(
  data.table::as.data.table(labels$cluster_map),
  file.path(result_dir, "cluster_annotations.tsv"),
  sep = "\t"
)
data.table::fwrite(
  data.table::as.data.table(cell_scored$coverage),
  file.path(result_dir, "cell_type_marker_coverage.tsv"),
  sep = "\t"
)
data.table::fwrite(
  data.table::as.data.table(state_scored$coverage),
  file.path(result_dir, "state_marker_coverage.tsv"),
  sep = "\t"
)

umap <- as.data.frame(SingleCellExperiment::reducedDim(sce, "UMAP"))
names(umap) <- c("UMAP1", "UMAP2")
umap$cell_id <- colnames(sce)
umap$sample_id <- sce$sample_id
umap$cluster_id <- as.character(sce$cluster_id)
umap$coarse_cell_type <- sce$coarse_cell_type
data.table::fwrite(
  data.table::as.data.table(umap),
  file.path(result_dir, "umap_coordinates.tsv.gz"),
  sep = "\t"
)

write_session_manifest(
  file.path(result_dir, "analysis_manifest.json"),
  list(
    dataset_id = dataset_id,
    input = normalizePath(input),
    output = normalizePath(output, mustWork = FALSE),
    retained_cells = ncol(sce),
    hvg_count = length(hvg),
    latent_method = paste(
      "Uncorrected PCA. Sample/batch was not regressed because capture",
      "library is confounded with biological time in core time courses."
    ),
    annotation_backend = paste(
      "Conservative cluster-level canonical marker scoring with explicit",
      "unknown/ambiguous class"
    ),
    state_guardrail = paste(
      "State scores are provisional expression evidence axes and are not",
      "interpreted as fate without temporal, lineage, spatial or perturbation",
      "support."
    )
  )
)
message("Completed annotation for ", dataset_id)

  },
  "integrate" = function() {

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
data_root <- Sys.getenv("MUSCLE_DATA_ROOT", "data")
processed_root <- Sys.getenv(
  "MUSCLE_PROCESSED_ROOT", file.path(data_root, "processed")
)
input_manifest <- if (length(args) >= 1L) args[[1]] else
  "config/core_state_space_inputs.tsv"
output <- if (length(args) >= 2L) args[[2]] else
  file.path(processed_root, "integrated/core_regeneration_state_space.rds")
reference_cells_per_dataset <- if (length(args) >= 3L) {
  as.integer(args[[3]])
} else {
  1200L
}

assert_packages(c(
  "SingleCellExperiment", "SummarizedExperiment", "Matrix", "matrixStats",
  "irlba", "harmony", "uwot", "igraph", "scran", "data.table",
  "S4Vectors", "jsonlite"
))
if (!file.exists(input_manifest)) stop("Missing input manifest: ", input_manifest)
manifest <- data.table::fread(input_manifest, sep = "\t")
manifest[, path := sub(
  "^\\$\\{MUSCLE_PROCESSED_ROOT\\}",
  processed_root, path
)]
manifest <- manifest[file.exists(path)]
if (nrow(manifest) < 3L) {
  stop("At least three completed annotated datasets are required")
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
result_dir <- "results/03_state_space"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(18018)

message("Auditing shared genes and collecting compact metadata")
gene_lists <- vector("list", nrow(manifest))
meta_list <- vector("list", nrow(manifest))
names(gene_lists) <- names(meta_list) <- manifest$dataset_id

for (i in seq_len(nrow(manifest))) {
  sce <- readRDS(manifest$path[[i]])
  gene_lists[[i]] <- rownames(sce)
  meta <- data.table::as.data.table(
    as.data.frame(SummarizedExperiment::colData(sce)),
    keep.rownames = "source_cell_id"
  )
  meta[, dataset_id := manifest$dataset_id[[i]]]
  meta[, injury_model := manifest$model[[i]]]
  meta[, integrated_cell_id := paste(dataset_id, source_cell_id, sep = "::")]
  meta_list[[i]] <- meta
  rm(sce)
  invisible(gc())
}
common_genes <- Reduce(intersect, gene_lists)
if (length(common_genes) < 10000L) {
  stop("Too few genes shared across datasets: ", length(common_genes))
}

balanced_reference <- function(meta, n_target) {
  by_sample <- split(seq_len(nrow(meta)), meta$sample_id)
  per_sample <- max(1L, ceiling(n_target / length(by_sample)))
  picked <- unlist(lapply(by_sample, function(ix) {
    sample(ix, min(length(ix), per_sample))
  }), use.names = FALSE)
  if (length(picked) > n_target) picked <- sample(picked, n_target)
  sort(picked)
}

message("Building a balanced cross-dataset reference matrix")
reference_matrices <- vector("list", nrow(manifest))
reference_index <- vector("list", nrow(manifest))
names(reference_matrices) <- names(reference_index) <- manifest$dataset_id
for (i in seq_len(nrow(manifest))) {
  sce <- readRDS(manifest$path[[i]])
  ix <- balanced_reference(meta_list[[i]], reference_cells_per_dataset)
  reference_index[[i]] <- ix
  reference_matrices[[i]] <- as.matrix(
    SummarizedExperiment::assay(sce, "logcounts")[
      common_genes, ix, drop = FALSE
    ]
  )
  rm(sce)
  invisible(gc())
}
reference <- do.call(cbind, reference_matrices)
rm(reference_matrices)
invisible(gc())

variance <- matrixStats::rowVars(reference)
nuisance <- grepl(
  "^(mt-|Rpl|Rps|Hba-|Hbb-|Gm[0-9]+)",
  common_genes,
  ignore.case = TRUE
)
hvg <- common_genes[order(variance, decreasing = TRUE)]
hvg <- head(hvg[!nuisance[match(hvg, common_genes)]], 2500L)
reference_hvg <- t(reference[hvg, , drop = FALSE])
rm(reference)
invisible(gc())

message("Learning shared uncorrected PCA on balanced cells")
pca <- irlba::prcomp_irlba(
  reference_hvg,
  n = 50,
  center = TRUE,
  scale. = TRUE
)
rm(reference_hvg)
invisible(gc())

project_chunks <- function(path, genes, pca_fit, chunk_size = 4000L) {
  sce <- readRDS(path)
  logcounts <- SummarizedExperiment::assay(sce, "logcounts")
  out <- matrix(
    NA_real_,
    nrow = ncol(sce),
    ncol = ncol(pca_fit$rotation),
    dimnames = list(colnames(sce), colnames(pca_fit$rotation))
  )
  starts <- seq.int(1L, ncol(sce), by = chunk_size)
  for (start in starts) {
    end <- min(ncol(sce), start + chunk_size - 1L)
    message("Projecting ", basename(path), " cells ", start, "-", end)
    x <- t(as.matrix(logcounts[genes, start:end, drop = FALSE]))
    x <- scale(x, center = pca_fit$center, scale = pca_fit$scale)
    x[!is.finite(x)] <- 0
    out[start:end, ] <- x %*% pca_fit$rotation
  }
  rm(sce)
  invisible(gc())
  out
}

pca_all <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  pca_all[[i]] <- project_chunks(manifest$path[[i]], hvg, pca)
}
pca_all <- do.call(rbind, pca_all)
metadata <- data.table::rbindlist(meta_list, fill = TRUE, use.names = TRUE)
if (nrow(metadata) != nrow(pca_all)) stop("Metadata/PCA row mismatch")
rownames(pca_all) <- metadata$integrated_cell_id

message("Aligning technical dataset effects with Harmony")
harmony_embedding <- harmony::RunHarmony(
  data_mat = pca_all,
  meta_data = as.data.frame(metadata[, .(dataset_id)]),
  vars_use = "dataset_id",
  max_iter = 30,
  return_object = FALSE,
  verbose = TRUE
)
rownames(harmony_embedding) <- metadata$integrated_cell_id

message("Constructing shared neighborhoods and UMAP")
graph <- scran::buildSNNGraph(harmony_embedding, k = 30, transposed = TRUE)
cluster <- igraph::cluster_louvain(graph)$membership
umap <- uwot::umap(
  harmony_embedding,
  n_neighbors = 30,
  min_dist = 0.25,
  metric = "cosine",
  n_threads = 1,
  verbose = TRUE,
  ret_model = FALSE
)
colnames(umap) <- c("UMAP1", "UMAP2")
rownames(umap) <- metadata$integrated_cell_id

metadata[, integrated_cluster := as.character(cluster)]
cluster_audit <- metadata[, {
  type_tab <- table(coarse_cell_type)
  dataset_tab <- table(dataset_id)
  p <- as.numeric(dataset_tab) / sum(dataset_tab)
  list(
    n_cells = .N,
    n_datasets = data.table::uniqueN(dataset_id),
    dataset_entropy = -sum(p * log(p + 1e-12)),
    normalized_dataset_entropy = -sum(p * log(p + 1e-12)) /
      log(length(dataset_tab)),
    min_dataset_fraction = min(p),
    max_dataset_fraction = max(p),
    majority_cell_type = names(type_tab)[which.max(type_tab)],
    cell_type_concordance = max(type_tab) / sum(type_tab)
  )
}, by = integrated_cluster]
metadata <- merge(
  metadata,
  cluster_audit[, .(
    integrated_cluster, integrated_majority_cell_type = majority_cell_type,
    integrated_cell_type_concordance = cell_type_concordance
  )],
  by = "integrated_cluster",
  all.x = TRUE,
  sort = FALSE
)
metadata <- metadata[match(rownames(pca_all), integrated_cell_id)]
metadata[, unified_cell_type := ifelse(
  integrated_cell_type_concordance >= 0.65,
  integrated_majority_cell_type,
  "Unknown_or_cross_dataset_discordant"
)]

# This deliberately remains an evidence label. Fate labels are introduced only
# after adjacent-time persistence, lineage, spatial and perturbation checks.
metadata[, unified_state_evidence := provisional_state_axis]
metadata[, fate_interpretation_status := ifelse(
  unified_state_evidence == "Activation",
  "transient_stress_or_activation_unresolved",
  "provisional_expression_axis_not_fate"
)]

compact_col_data <- S4Vectors::DataFrame(as.data.frame(metadata))
rownames(compact_col_data) <- metadata$integrated_cell_id
compact <- SingleCellExperiment::SingleCellExperiment(colData = compact_col_data)
SingleCellExperiment::reducedDim(compact, "PCA_UNCORRECTED") <- pca_all
SingleCellExperiment::reducedDim(compact, "HARMONY") <- harmony_embedding
SingleCellExperiment::reducedDim(compact, "UMAP") <- umap
base::saveRDS(compact, output, compress = FALSE)

state_score_cols <- grep("^state_score_", names(metadata), value = TRUE)
keep_cols <- c(
  "integrated_cell_id", "dataset_id", "injury_model", "sample_id",
  "time_days", "coarse_cell_type", "unified_cell_type",
  "unified_state_evidence", "fate_interpretation_status",
  "integrated_cluster", state_score_cols
)
optional_design_cols <- intersect(
  c("age_group", "condition", "sex", "age", "genotype", "replicate"),
  names(metadata)
)
keep_cols <- c(keep_cols, optional_design_cols)
umap_dt <- metadata[, ..keep_cols]
umap_dt[, `:=`(UMAP1 = umap[, 1], UMAP2 = umap[, 2])]
data.table::fwrite(
  umap_dt, file.path(result_dir, "integrated_cell_state_metadata.tsv.gz"),
  sep = "\t"
)
data.table::fwrite(
  cluster_audit, file.path(result_dir, "integrated_cluster_audit.tsv"),
  sep = "\t"
)
data.table::fwrite(
  data.table::data.table(gene = hvg, rank = seq_along(hvg)),
  file.path(result_dir, "shared_hvg.tsv"),
  sep = "\t"
)
data.table::fwrite(
  manifest, file.path(result_dir, "included_datasets.tsv"), sep = "\t"
)

write_session_manifest(
  file.path(result_dir, "analysis_manifest.json"),
  list(
    input_manifest = normalizePath(input_manifest),
    output = normalizePath(output, mustWork = FALSE),
    datasets = manifest$dataset_id,
    n_cells = ncol(compact),
    common_genes = length(common_genes),
    hvg_count = length(hvg),
    reference_cells_per_dataset = reference_cells_per_dataset,
    uncorrected_embedding = paste(
      "Shared PCA learned from a dataset-balanced reference; retained for",
      "leave-dataset and overcorrection sensitivity analyses."
    ),
    aligned_embedding = paste(
      "Harmony on dataset identity only; biological time, model, age and",
      "condition were not regressed."
    ),
    state_guardrail = paste(
      "Unified states are expression evidence axes. Activation is explicitly",
      "unresolved between transient stress and biological activation; fate",
      "requires temporal persistence plus lineage/spatial/perturbation evidence."
    )
  )
)
message("Completed core Muscle Regeneration State Space")

  },
  "discovery" = function() {
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(data.table)
  library(ggplot2)
  library(mgcv)
})

input <- file.path(Sys.getenv("MUSCLE_PROCESSED_ROOT", "data/processed"), "integrated/core_regeneration_state_space.rds")
outdir <- "results/discovery_clock"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input))

weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

soft_phase <- function(m, centers, temperature = 0.8) {
  z <- sweep(m, 1L, apply(m, 1L, max, na.rm = TRUE), "-") / temperature
  p <- exp(z)
  p <- p / rowSums(p)
  as.numeric(p %*% centers)
}

sce <- readRDS(input)
d <- as.data.table(as.data.frame(colData(sce)))
rm(sce); gc()

required <- c("dataset_id", "sample_id", "time_days", "passes_qc_and_singlet",
              "unified_cell_type", "state_score_Activation", "state_score_Proliferation",
              "state_score_Commitment", "state_score_Fusion", "state_score_Remodeling",
              "state_score_Resolution", "age_group", "replicate", "condition")
stopifnot(all(required %in% names(d)))
d <- d[passes_qc_and_singlet %in% TRUE]
d <- d[!grepl("^Unknown", unified_cell_type)]

type_map <- c(
  Myogenic_progenitor = "MuSC_myogenic", Myonucleus = "myonucleus",
  FAP = "FAP", Macrophage = "macrophage", Neutrophil = "neutrophil",
  Endothelial = "endothelial", Pericyte = "pericyte", T_cell = "T_NK",
  B_cell = "other_immune", Schwann = "other_major"
)
d[, cell_type := unname(type_map[unified_cell_type])]
d <- d[!is.na(cell_type)]

score_cols <- c("state_score_Activation", "state_score_Proliferation",
                "state_score_Commitment", "state_score_Fusion",
                "state_score_Remodeling", "state_score_Resolution")
centers <- c(0.12, 0.28, 0.45, 0.64, 0.80, 0.94)
d[, provisional_phase := soft_phase(as.matrix(.SD), centers), .SDcols = score_cols]

d[, age_class := fifelse(dataset_id == "GSE232106", as.character(age_group),
                  fifelse(dataset_id %in% c("GSE138826", "GSE143437"), "young", NA_character_))]
d[, injury_model_std := fifelse(dataset_id == "GSE138826", "CTX",
                         fifelse(dataset_id %in% c("GSE143437", "GSE232106"), "NTX",
                         fifelse(dataset_id == "GSE215808", "VML", "unknown")))]

# A/B and C/D are paired-TA capture groups for geriatric GSE232106. This is a
# conservative inferred grouping and is kept explicit in the output manifest.
d[, biological_group := sample_id]
d[dataset_id == "GSE232106" & age_class == "geriatric",
  biological_group := paste(dataset_id, age_class, time_days,
                            fifelse(replicate %in% c("A", "B"), "mouse_pair_1",
                            fifelse(replicate %in% c("C", "D"), "mouse_pair_2", replicate)), sep = "::")]

cell_summary <- d[, c(
  list(n_cells = .N,
       phase = median(provisional_phase, na.rm = TRUE),
       phase_se = 1.2533 * mad(provisional_phase, na.rm = TRUE) / sqrt(.N)),
  lapply(.SD, median, na.rm = TRUE)),
  by = .(dataset_id, sample_id, biological_group, time_days, age_class,
         injury_model_std, condition, cell_type), .SDcols = score_cols]
cell_summary[, reliability := pmin(1, sqrt(n_cells / 100)) / pmax(phase_se, 0.02)]
cell_summary[, reliability := reliability / max(reliability, na.rm = TRUE), by = cell_type]

sample_clock <- cell_summary[, .(
  global_tau = weighted_median(phase, reliability),
  effective_cell_types = sum(n_cells >= 20),
  total_cells = sum(n_cells)),
  by = .(dataset_id, sample_id, biological_group, time_days, age_class,
         injury_model_std, condition)]
cell_summary <- merge(cell_summary, sample_clock,
                      by = c("dataset_id", "sample_id", "biological_group", "time_days",
                             "age_class", "injury_model_std", "condition"))

reference <- cell_summary[age_class == "young" & injury_model_std == "NTX" & n_cells >= 20]
stopifnot(uniqueN(reference$time_days) >= 4L)
grid <- seq(0, 1, length.out = 201)

# Expected normal lineage phase at consensus biological time. Low-cell samples
# are excluded from fitting, and each cell type is fitted separately.
reference_models <- list()
for (ct in unique(reference$cell_type)) {
  z <- reference[cell_type == ct]
  if (nrow(z) < 6L || uniqueN(z$global_tau) < 4L) next
  k <- min(5L, uniqueN(z$global_tau) - 1L)
  reference_models[[ct]] <- gam(phase ~ s(global_tau, k = k), data = z,
                                 weights = reliability, method = "REML")
}
cell_summary[, expected_phase := NA_real_]
for (ct in names(reference_models)) {
  idx <- which(cell_summary$cell_type == ct)
  cell_summary$expected_phase[idx] <- predict(reference_models[[ct]],
                                               newdata = cell_summary[idx])
}
cell_summary[, deviation := phase - expected_phase]

ref_dev <- cell_summary[age_class == "young" & injury_model_std == "NTX" &
                          n_cells >= 20 & is.finite(deviation), deviation]
sigma_reference <- max(0.05, sqrt(mean(ref_dev^2, na.rm = TRUE)))

rsi <- cell_summary[n_cells >= 20 & is.finite(deviation), .(
  RSI = regenerative_synchrony_index(deviation, reliability,
                                      sigma_reference = sigma_reference),
  phase_dispersion = weighted_var(deviation, reliability),
  effective_cell_types = .N),
  by = .(dataset_id, sample_id, biological_group, time_days, age_class,
         injury_model_std, condition)]
rsi[, reference_sigma := sigma_reference]

# Aggregate capture libraries to the conservative biological group before age contrasts.
bio_rsi <- rsi[, .(RSI = weighted.mean(RSI, effective_cell_types, na.rm = TRUE),
                   phase_dispersion = weighted.mean(phase_dispersion, effective_cell_types,
                                                    na.rm = TRUE),
                   n_capture_libraries = .N),
               by = .(dataset_id, biological_group, time_days, age_class,
                      injury_model_std, condition)]

age_overlap <- bio_rsi[dataset_id == "GSE232106" & time_days %in% c(1, 3.5) &
                         age_class %in% c("young", "geriatric")]
age_effects <- age_overlap[, .(
  n_units = .N, mean_RSI = mean(RSI), sd_RSI = sd(RSI)),
  by = .(time_days, age_class)]
age_wide <- dcast(age_effects, time_days ~ age_class, value.var = "mean_RSI")
if (all(c("young", "geriatric") %in% names(age_wide)))
  age_wide[, geriatric_minus_young := geriatric - young]

# Prespecified directional phase lags.
lag_pairs <- data.table(a = c("MuSC_myogenic", "FAP", "endothelial"),
                        b = c("macrophage", "MuSC_myogenic", "myonucleus"),
                        contrast = c("MuSC-minus-macrophage", "FAP-minus-MuSC",
                                     "endothelial-minus-myonucleus"))
phase_wide <- dcast(cell_summary[n_cells >= 20],
                    dataset_id + sample_id + biological_group + time_days + age_class +
                      injury_model_std + condition ~ cell_type,
                    value.var = "phase")
lags <- rbindlist(lapply(seq_len(nrow(lag_pairs)), function(i) {
  p <- lag_pairs[i]
  if (!all(c(p$a, p$b) %in% names(phase_wide))) return(NULL)
  phase_wide[, .(dataset_id, sample_id, biological_group, time_days, age_class,
                 injury_model_std, condition, contrast = p$contrast,
                 phase_lag = get(p$a) - get(p$b))]
}))

# Exploratory open-phase coupling summary. Ridge regularization stabilizes the
# underdetermined short time series; coefficients are predictive, not causal.
coupling_types <- c("MuSC_myogenic", "macrophage", "FAP", "endothelial", "myonucleus")
young_curve <- cell_summary[age_class == "young" & injury_model_std == "NTX" &
                              cell_type %in% coupling_types & n_cells >= 20,
                            .(phase = weighted.mean(phase, reliability)),
                            by = .(time_days, cell_type)]
curve_wide <- dcast(young_curve, time_days ~ cell_type, value.var = "phase")
curve_wide <- curve_wide[complete.cases(curve_wide[, ..coupling_types])][order(time_days)]
K <- matrix(NA_real_, length(coupling_types), length(coupling_types),
            dimnames = list(target = coupling_types, source = coupling_types))
if (nrow(curve_wide) >= 5L) {
  theta <- as.matrix(curve_wide[, ..coupling_types]) * pi
  dt <- diff(curve_wide$time_days)
  deriv <- sweep(diff(theta), 1L, dt, "/")
  mid <- (theta[-1, , drop = FALSE] + theta[-nrow(theta), , drop = FALSE]) / 2
  lambda <- 0.25
  for (i in seq_along(coupling_types)) {
    others <- setdiff(seq_along(coupling_types), i)
    X <- cbind(omega = 1, sapply(others, function(j) sin(mid[, j] - mid[, i])))
    pen <- diag(c(0, rep(lambda, length(others))))
    beta <- solve(crossprod(X) + pen, crossprod(X, deriv[, i]))
    K[i, others] <- beta[-1]
    K[i, i] <- 0
  }
}

fwrite(cell_summary, file.path(outdir, "sample_celltype_phase.tsv"), sep = "\t")
fwrite(rsi, file.path(outdir, "sample_RSI.tsv"), sep = "\t")
fwrite(bio_rsi, file.path(outdir, "biological_unit_RSI.tsv"), sep = "\t")
fwrite(age_effects, file.path(outdir, "age_overlap_summary.tsv"), sep = "\t")
if (exists("age_wide")) fwrite(age_wide, file.path(outdir, "age_overlap_effect.tsv"), sep = "\t")
fwrite(lags, file.path(outdir, "prespecified_phase_lags.tsv"), sep = "\t")
fwrite(as.data.table(as.table(K)), file.path(outdir, "exploratory_coupling_K.tsv"), sep = "\t")
writeLines(c(
  paste("input:", input),
  paste("cells_after_existing_QC_and_annotation_filter:", nrow(d)),
  paste("samples:", nrow(rsi)),
  paste("biological_groups:", nrow(bio_rsi)),
  paste("sigma_reference:", signif(sigma_reference, 5)),
  "phase_estimator: softmax-weighted ordered expression-state programs",
  "inference_unit: biological_group; cells are measurement units only",
  "coupling_status: exploratory ridge-stabilized open-phase summary; not causal"
), file.path(outdir, "analysis_manifest.txt"))

p1 <- ggplot(bio_rsi[!is.na(age_class) & is.finite(RSI)],
             aes(time_days, RSI, colour = age_class)) +
  geom_point(size = 1.7) +
  facet_wrap(~dataset_id, scales = "free_x") +
  scale_colour_manual(values = c(young = "#2878B5", geriatric = "#D1495B")) +
  labs(x = "Days post injury", y = "Regenerative Synchrony Index", colour = "Age") +
  theme_classic(base_size = 9)
ggsave(file.path(outdir, "RSI_trajectories.pdf"), p1, width = 7.2, height = 3.7)

p2 <- ggplot(cell_summary[dataset_id == "GSE232106" & time_days %in% c(1, 3.5) &
                            age_class %in% c("young", "geriatric") & n_cells >= 20 &
                            is.finite(deviation)],
             aes(cell_type, deviation, colour = age_class)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
             size = 0.8, alpha = 0.65) +
  facet_wrap(~time_days, labeller = label_both) +
  scale_colour_manual(values = c(young = "#2878B5", geriatric = "#D1495B")) +
  labs(x = NULL, y = "Phase deviation from young reference", colour = "Age") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(outdir, "age_phase_deviations.pdf"), p2, width = 8, height = 4.2)

cat("Discovery clock complete:", nrow(d), "cells,", nrow(rsi), "samples,",
    nrow(bio_rsi), "biological groups.\n")

  },
  "multiview" = function() {
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

indir <- "results/discovery_clock"
outdir <- "results/multiview_clock"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(indir, "sample_celltype_phase.tsv"))

score_cols <- c("state_score_Activation", "state_score_Proliferation",
                "state_score_Commitment", "state_score_Fusion",
                "state_score_Remodeling", "state_score_Resolution")
stopifnot(all(score_cols %in% names(x)))
x[, target_tau := pmin(1, log1p(time_days) / log(8))]
x[, row_id := .I]
reference <- x[age_class == "young" & injury_model_std %in% c("NTX", "CTX") &
                 n_cells >= 20 & complete.cases(x[, ..score_cols])]

ridge_fit <- function(train, features, lambda = 1) {
  X <- as.matrix(train[, ..features])
  mu <- colMeans(X); sc <- apply(X, 2, sd)
  sc[!is.finite(sc) | sc < 1e-8] <- 1
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sc, "/")
  D <- diag(c(0, rep(lambda, ncol(Z))))
  gram <- crossprod(cbind(1, Z))
  if (!identical(dim(gram), dim(D))) {
    stop(sprintf("ridge dimension mismatch: n=%d p=%d gram=%s penalty=%s features=%s",
                 nrow(Z), ncol(Z), paste(dim(gram), collapse = "x"),
                 paste(dim(D), collapse = "x"), paste(features, collapse = ",")))
  }
  beta <- solve(gram + D,
                crossprod(cbind(1, Z), train$target_tau))
  list(beta = beta, mu = mu, sc = sc, features = features)
}

ridge_predict <- function(fit, test) {
  features <- fit$features
  X <- as.matrix(test[, ..features])
  Z <- sweep(sweep(X, 2, fit$mu, "-"), 2, fit$sc, "/")
  pmin(1.25, pmax(-0.15, as.numeric(cbind(1, Z) %*% fit$beta)))
}

choose_lambda <- function(train, features, groups = "sample_id",
                          grid = c(0.01, 0.1, 1, 10, 100)) {
  ids <- unique(train[[groups]])
  if (length(ids) < 5L) return(1)
  err <- sapply(grid, function(lam) {
    pred <- rep(NA_real_, nrow(train))
    for (id in ids) {
      te <- which(train[[groups]] == id); tr <- setdiff(seq_len(nrow(train)), te)
      if (length(tr) <= length(features) + 1L) next
      fit <- ridge_fit(train[tr], features, lam)
      pred[te] <- ridge_predict(fit, train[te])
    }
    mean(abs(pred - train$target_tau), na.rm = TRUE)
  })
  if (!any(is.finite(err))) return(1)
  err[!is.finite(err)] <- Inf
  grid[which.min(err)]
}

weighted_median <- function(v, w) {
  ok <- is.finite(v) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  o <- order(v[ok]); vv <- v[ok][o]; ww <- w[ok][o]
  vv[which(cumsum(ww) >= sum(ww) / 2)[1]]
}

# Nested leave-study-out comparison. Every fitted coefficient and calibration
# is learned without the held-out study.
studies <- unique(reference$dataset_id)
cv <- list()
for (held in studies) {
  train <- reference[dataset_id != held]
  test <- reference[dataset_id == held]
  if (nrow(train) < 15L || nrow(test) == 0L) next

  null_pred <- rep(mean(train$target_tau), nrow(test))
  lam_pool <- choose_lambda(train, score_cols)
  pooled_pred <- ridge_predict(ridge_fit(train, score_cols, lam_pool), test)

  pseudo_feature <- "phase"
  lam_pseudo <- choose_lambda(train, pseudo_feature)
  pseudo_pred <- ridge_predict(ridge_fit(train, pseudo_feature, lam_pseudo), test)

  ind_pred <- rep(NA_real_, nrow(test))
  for (ct in unique(test$cell_type)) {
    tr <- train[cell_type == ct]; te <- which(test$cell_type == ct)
    if (nrow(tr) < 8L || uniqueN(tr$time_days) < 4L) next
    lam <- choose_lambda(tr, score_cols)
    ind_pred[te] <- ridge_predict(ridge_fit(tr, score_cols, lam), test[te])
  }
  tmp <- copy(test)
  tmp[, independent_phase := ind_pred]
  coupled <- tmp[, .(pred = weighted_median(independent_phase, reliability),
                     target_tau = unique(target_tau)[1]), by = sample_id]

  sample_base <- unique(test[, .(sample_id, target_tau)])
  pool_s <- data.table(sample_id = test$sample_id, p = pooled_pred)[,
             .(pred = median(p)), by = sample_id]
  pseudo_s <- data.table(sample_id = test$sample_id, p = pseudo_pred)[,
               .(pred = median(p)), by = sample_id]
  null_s <- data.table(sample_id = test$sample_id, p = null_pred)[,
             .(pred = median(p)), by = sample_id]
  ind_s <- tmp[, .(pred = weighted.mean(independent_phase, reliability, na.rm = TRUE)),
               by = sample_id]

  add_metric <- function(tab, model) {
    z <- merge(sample_base, tab, by = "sample_id", all.x = TRUE)
    z[, .(held_out_study = held, model = model, sample_id,
          target_tau, prediction = pred, absolute_error = abs(pred - target_tau))]
  }
  cv[[paste(held, "null")]] <- add_metric(null_s, "absolute_time_null")
  cv[[paste(held, "pooled")]] <- add_metric(pool_s, "cell_state_pooled")
  cv[[paste(held, "pseudo")]] <- add_metric(pseudo_s, "global_pseudotime")
  cv[[paste(held, "ind")]] <- add_metric(ind_s, "independent_celltype_clocks")
  cv[[paste(held, "coupled")]] <- add_metric(coupled[, .(sample_id, pred)],
                                              "coupled_regeneration_clock")
}
cv <- rbindlist(cv, fill = TRUE)
cv_summary <- cv[, .(n_samples = sum(is.finite(absolute_error)),
                     MAE_tau = mean(absolute_error, na.rm = TRUE),
                     median_AE_tau = median(absolute_error, na.rm = TRUE)),
                 by = .(held_out_study, model)]
cv_overall <- cv[, .(n_samples = sum(is.finite(absolute_error)),
                     MAE_tau = mean(absolute_error, na.rm = TRUE),
                     median_AE_tau = median(absolute_error, na.rm = TRUE)), by = model]

# Final clocks for biological contrasts; trained only in young successful repair.
x[, phase_mv := NA_real_]
lambda_by_type <- list()
for (ct in unique(reference$cell_type)) {
  tr <- reference[cell_type == ct]
  if (nrow(tr) < 8L || uniqueN(tr$time_days) < 4L) next
  lam <- choose_lambda(tr, score_cols)
  lambda_by_type[[ct]] <- lam
  idx <- which(x$cell_type == ct & complete.cases(x[, ..score_cols]))
  x$phase_mv[idx] <- ridge_predict(ridge_fit(tr, score_cols, lam), x[idx])
}

sample_tau <- x[is.finite(phase_mv) & n_cells >= 20,
                .(consensus_tau = weighted_median(phase_mv, reliability),
                  effective_cell_types = .N),
                by = .(dataset_id, sample_id, biological_group, time_days, age_class,
                       injury_model_std, condition)]
x <- merge(x, sample_tau, by = c("dataset_id", "sample_id", "biological_group",
                                  "time_days", "age_class", "injury_model_std", "condition"),
           all.x = TRUE)

# Learn normal, nonzero lineage offsets as a function of consensus time.
x[, expected_mv := NA_real_]
ref2 <- x[age_class == "young" & injury_model_std %in% c("NTX", "CTX") &
            is.finite(phase_mv) & is.finite(consensus_tau) & n_cells >= 20]
for (ct in unique(ref2$cell_type)) {
  tr <- ref2[cell_type == ct]
  if (nrow(tr) < 8L || uniqueN(tr$consensus_tau) < 4L) next
  fit <- lm(phase_mv ~ poly(consensus_tau, 2, raw = TRUE), data = tr,
            weights = reliability)
  idx <- which(x$cell_type == ct & is.finite(x$consensus_tau))
  x$expected_mv[idx] <- predict(fit, newdata = x[idx])
}
x[, deviation_mv := phase_mv - expected_mv]
sigma_ref <- max(0.05, sqrt(mean(ref2$phase_mv^2 * 0 +
  x[match(ref2$row_id, x$row_id), deviation_mv]^2, na.rm = TRUE)))

rsi <- x[is.finite(deviation_mv) & n_cells >= 20,
         .(RSI_mv = regenerative_synchrony_index(deviation_mv, reliability,
                                                  sigma_reference = sigma_ref),
           phase_dispersion_mv = weighted_var(deviation_mv, reliability),
           effective_cell_types = .N),
         by = .(dataset_id, sample_id, biological_group, time_days, age_class,
                injury_model_std, condition)]
bio <- rsi[, .(RSI_mv = weighted.mean(RSI_mv, effective_cell_types, na.rm = TRUE),
               phase_dispersion_mv = weighted.mean(phase_dispersion_mv,
                                                   effective_cell_types, na.rm = TRUE),
               n_capture_libraries = .N),
           by = .(dataset_id, biological_group, time_days, age_class,
                  injury_model_std, condition)]

overlap <- bio[dataset_id == "GSE232106" & time_days %in% c(1, 3.5) &
                 age_class %in% c("young", "geriatric")]
age_summary <- overlap[, .(n_units = .N, mean_RSI = mean(RSI_mv),
                           sd_RSI = sd(RSI_mv)), by = .(time_days, age_class)]
age_effect <- dcast(age_summary, time_days ~ age_class, value.var = "mean_RSI")
if (all(c("young", "geriatric") %in% names(age_effect)))
  age_effect[, geriatric_minus_young := geriatric - young]

fwrite(x, file.path(outdir, "multiview_celltype_clock.tsv"), sep = "\t")
fwrite(rsi, file.path(outdir, "sample_RSI_multiview.tsv"), sep = "\t")
fwrite(bio, file.path(outdir, "biological_unit_RSI_multiview.tsv"), sep = "\t")
fwrite(cv, file.path(outdir, "leave_study_out_predictions.tsv"), sep = "\t")
fwrite(cv_summary, file.path(outdir, "leave_study_out_metrics.tsv"), sep = "\t")
fwrite(cv_overall, file.path(outdir, "model_comparison_overall.tsv"), sep = "\t")
fwrite(age_summary, file.path(outdir, "age_overlap_summary.tsv"), sep = "\t")
fwrite(age_effect, file.path(outdir, "age_overlap_effect.tsv"), sep = "\t")
fwrite(data.table(cell_type = names(lambda_by_type),
                  lambda = unlist(lambda_by_type)),
       file.path(outdir, "ridge_lambda_by_celltype.tsv"), sep = "\t")

p <- ggplot(cv_overall, aes(reorder(model, MAE_tau), MAE_tau)) +
  geom_col(fill = "#4477AA") + coord_flip() +
  labs(x = NULL, y = "Leave-study-out biological-time MAE") +
  theme_classic(base_size = 9)
ggsave(file.path(outdir, "model_comparison.pdf"), p, width = 5.3, height = 3.3)

p2 <- ggplot(bio[!is.na(age_class) & is.finite(RSI_mv)],
             aes(time_days, RSI_mv, colour = age_class)) +
  geom_point(position = position_jitter(width = 0.08), size = 1.6) +
  facet_wrap(~dataset_id, scales = "free_x") +
  scale_colour_manual(values = c(young = "#2878B5", geriatric = "#D1495B")) +
  labs(x = "Days post injury", y = "Multiview RSI", colour = "Age") +
  theme_classic(base_size = 9)
ggsave(file.path(outdir, "multiview_RSI.pdf"), p2, width = 7.2, height = 3.7)

writeLines(c(
  paste("young_reference_rows:", nrow(reference)),
  paste("reference_studies:", paste(studies, collapse = ",")),
  paste("sigma_reference:", signif(sigma_ref, 5)),
  "validation: coefficients and calibration refit inside each held-out-study fold",
  "boundary: reference views are expression-state modules only; regulon/pathway views await raw-count scoring"
), file.path(outdir, "analysis_manifest.txt"))

cat("Multiview clock and nested leave-study-out comparison complete.\n")

  },
  "depletion" = function() {
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(scran)
  library(scater)
  library(bluster)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE291782")
object_out <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "processed/GSE291782_qc_annotated_sce.rds")
outdir <- "results/GSE291782_failure"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_out), recursive = TRUE, showWarnings = FALSE)

matrix_files <- list.files(input_dir, pattern = "^[^._].*_matrix\\.mtx\\.gz$", full.names = TRUE)
stopifnot(length(matrix_files) == 6L)

read_bundle <- function(mtx_file) {
  prefix <- sub("_matrix\\.mtx\\.gz$", "", mtx_file)
  feat <- fread(paste0(prefix, "_features.tsv.gz"), header = FALSE)
  bc <- fread(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)[[1]]
  mat <- readMM(gzfile(mtx_file))
  mat <- as(mat, "dgCMatrix")
  rownames(mat) <- make.unique(feat[[2]])
  sample <- basename(prefix)
  colnames(mat) <- paste(sample, bc, sep = "::")
  condition <- if (grepl("_CLL_", sample)) "clodronate" else "sham"
  replicate <- sub(".*_(CLL|SHAM)_", "", sample)
  sce <- SingleCellExperiment(list(counts = mat))
  colData(sce)$sample_id <- sample
  colData(sce)$condition <- condition
  colData(sce)$replicate <- replicate
  sce
}

set.seed(20260824)
objects <- lapply(matrix_files, read_bundle)
thresholds <- list(); qc_rows <- list()

for (i in seq_along(objects)) {
  sce <- objects[[i]]
  mt <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
  qc <- perCellQCMetrics(sce, subsets = list(Mt = mt))
  sid <- unique(colData(sce)$sample_id)
  low_detected <- max(200, median(qc$detected) - 3 * mad(qc$detected))
  low_sum <- max(500, median(qc$sum) - 3 * mad(qc$sum))
  high_detected <- median(qc$detected) + 5 * mad(qc$detected)
  high_sum <- median(qc$sum) + 5 * mad(qc$sum)
  high_mt <- min(40, median(qc$subsets_Mt_percent) + 3 * mad(qc$subsets_Mt_percent))
  pass <- qc$detected >= low_detected & qc$sum >= low_sum &
    qc$detected <= high_detected & qc$sum <= high_sum &
    qc$subsets_Mt_percent <= high_mt
  colData(sce) <- cbind(colData(sce), qc)
  colData(sce)$passes_core_qc <- pass

  dbl_class <- rep(NA_character_, ncol(sce)); dbl_score <- rep(NA_real_, ncol(sce))
  if (sum(pass) >= 100L) {
    called <- scDblFinder(sce[, pass], samples = rep(sid, sum(pass)),
                          BPPARAM = SerialParam(), verbose = FALSE)
    dbl_class[pass] <- as.character(called$scDblFinder.class)
    dbl_score[pass] <- called$scDblFinder.score
  }
  colData(sce)$scDblFinder.class <- dbl_class
  colData(sce)$scDblFinder.score <- dbl_score
  colData(sce)$passes_qc_and_singlet <- pass & dbl_class == "singlet"
  objects[[i]] <- sce
  thresholds[[i]] <- data.table(sample_id = sid, low_detected, low_sum,
                                 high_detected, high_sum, high_mt)
  qc_rows[[i]] <- data.table(sample_id = sid, condition = unique(sce$condition),
                             input_cells = ncol(sce), core_qc_pass = sum(pass),
                             doublets = sum(dbl_class == "doublet", na.rm = TRUE),
                             final_singlets = sum(sce$passes_qc_and_singlet, na.rm = TRUE))
}

sce <- do.call(cbind, lapply(objects, function(z) z[, z$passes_qc_and_singlet %in% TRUE]))
rm(objects); gc()
sizeFactors(sce) <- NULL
sce <- logNormCounts(sce)
cat("Combined post-QC object:", nrow(sce), "genes x", ncol(sce), "cells; finite size factors:",
    sum(is.finite(sizeFactors(sce))), "\n")

# HVGs are used only to obtain a visualization/clustering space. Direct sparse
# row variances avoid modelGeneVar's residual-d.f. failure on this object.
gene_var <- DelayedMatrixStats::rowVars(logcounts(sce), na.rm = TRUE)
cat("Finite gene variances:", sum(is.finite(gene_var)), "of", length(gene_var), "\n")
hvg <- head(order(gene_var, decreasing = TRUE, na.last = NA),
            min(3000L, sum(is.finite(gene_var))))
if (length(hvg) < 50L) stop("HVG selection returned fewer than 50 finite-variance genes")
sce <- runPCA(sce, subset_row = hvg,
              ncomponents = min(30L, length(hvg) - 1L),
              BSPARAM = BiocSingular::IrlbaParam())
sce$cluster_id <- as.character(clusterRows(reducedDim(sce, "PCA"),
  NNGraphParam(k = 20, cluster.fun = "louvain")))

markers <- list(
  MuSC_myogenic = c("Pax7", "Calcr", "Sdc4", "Vcam1", "Myod1", "Myog", "Des"),
  myonucleus = c("Ttn", "Myh1", "Myh2", "Myh4", "Ckm", "Acta1"),
  FAP = c("Pdgfra", "Col1a1", "Col3a1", "Dcn", "Lum"),
  macrophage = c("Ptprc", "Lyz2", "Adgre1", "Cd68", "C1qa", "C1qb"),
  neutrophil = c("S100a8", "S100a9", "Retnlg", "Camp", "Ly6g"),
  T_NK = c("Cd3d", "Cd3e", "Trbc1", "Nkg7", "Klrd1"),
  endothelial = c("Pecam1", "Kdr", "Emcn", "Cdh5", "Klf2"),
  pericyte = c("Rgs5", "Cspg4", "Pdgfrb", "Notch3", "Des")
)
programs <- list(
  Activation = c("Fos", "Jun", "Junb", "Egr1", "Myod1"),
  Proliferation = c("Mki67", "Top2a", "Ccnb1", "Pcna", "Ube2c"),
  Commitment = c("Myod1", "Myf5", "Des", "Cdkn1c"),
  Fusion = c("Myog", "Mymk", "Mymx", "Tnnt1", "Ckm"),
  Remodeling = c("Col1a1", "Col3a1", "Fn1", "Mmp14", "Tgfb1"),
  Resolution = c("Mrc1", "C1qa", "C1qb", "Il10", "Apoe")
)

module_score <- function(genes) {
  hit <- intersect(genes, rownames(sce))
  if (length(hit) == 0L) return(rep(NA_real_, ncol(sce)))
  Matrix::colMeans(logcounts(sce)[hit, , drop = FALSE])
}
for (nm in names(markers)) colData(sce)[[paste0("marker_", nm)]] <- module_score(markers[[nm]])
for (nm in names(programs)) colData(sce)[[paste0("state_score_", nm)]] <- module_score(programs[[nm]])

md_df <- as.data.frame(colData(sce))
md <- as.data.table(md_df, keep.rownames = "cell_id")
marker_cols <- paste0("marker_", names(markers))
cluster_scores <- md[, lapply(.SD, mean, na.rm = TRUE), by = cluster_id, .SDcols = marker_cols]
M <- as.matrix(cluster_scores[, ..marker_cols])
top <- max.col(M, ties.method = "first")
ordered <- t(apply(M, 1, sort, decreasing = TRUE))
cluster_scores[, coarse_cell_type := names(markers)[top]]
cluster_scores[, annotation_margin := ordered[, 1] - ordered[, 2]]
cluster_scores[!is.finite(annotation_margin) | annotation_margin < 0.08,
               coarse_cell_type := "Unknown_or_ambiguous"]
md <- merge(md, cluster_scores[, .(cluster_id, coarse_cell_type, annotation_margin)],
            by = "cluster_id", sort = FALSE)

program_cols <- paste0("state_score_", names(programs))
sample_phase <- md[coarse_cell_type != "Unknown_or_ambiguous",
  c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
  by = .(sample_id, condition, replicate, coarse_cell_type), .SDcols = program_cols]
for (pc in program_cols) {
  mu <- mean(sample_phase[[pc]], na.rm = TRUE); ss <- sd(sample_phase[[pc]], na.rm = TRUE)
  sample_phase[[pc]] <- (sample_phase[[pc]] - mu) / ifelse(ss > 0, ss, 1)
}
centers <- c(0.12, 0.28, 0.45, 0.64, 0.80, 0.94)
soft_phase <- function(m) {
  z <- sweep(m, 1, apply(m, 1, max), "-") / 0.8
  p <- exp(z); p <- p / rowSums(p)
  as.numeric(p %*% centers)
}
sample_phase[, phase := soft_phase(as.matrix(.SD)), .SDcols = program_cols]

sample_phase[, expected_phase := NA_real_]
for (ct in unique(sample_phase$coarse_cell_type)) {
  sham_ids <- sample_phase[coarse_cell_type == ct & condition == "sham", sample_id]
  for (sid in sample_phase[coarse_cell_type == ct, unique(sample_id)]) {
    ref_ids <- if (sid %in% sham_ids) setdiff(sham_ids, sid) else sham_ids
    sample_phase[coarse_cell_type == ct & sample_id == sid,
                 expected_phase := mean(sample_phase[coarse_cell_type == ct &
                                                        sample_id %in% ref_ids, phase],
                                        na.rm = TRUE)]
  }
}
sample_phase[, deviation := phase - expected_phase]
sample_phase[, reliability := pmin(1, sqrt(n_cells / 100))]
sigma_ref <- max(0.05, sqrt(mean(sample_phase[condition == "sham", deviation]^2, na.rm = TRUE)))
rsi <- sample_phase[n_cells >= 20 & is.finite(deviation), .(
  RSI = regenerative_synchrony_index(deviation, reliability,
                                      sigma_reference = sigma_ref),
  dispersion = weighted_var(deviation, reliability),
  effective_cell_types = .N), by = .(sample_id, condition, replicate)]

rsi_valid <- rsi[is.finite(RSI)]
obs <- rsi_valid[, mean(RSI), by = condition]
observed_effect <- obs[condition == "clodronate", V1] - obs[condition == "sham", V1]
n_cll <- sum(rsi_valid$condition == "clodronate")
all_cll <- combn(rsi_valid$sample_id, n_cll, simplify = FALSE)
perm_effect <- vapply(all_cll, function(ids) {
  mean(rsi_valid[sample_id %in% ids, RSI]) - mean(rsi_valid[!sample_id %in% ids, RSI])
}, numeric(1))
perm_p_two <- mean(abs(perm_effect) >= abs(observed_effect) - 1e-12)
effect <- data.table(contrast = "clodronate_minus_sham", effect_RSI = observed_effect,
                     exact_permutation_p_two_sided = perm_p_two,
                     n_clodronate = sum(rsi_valid$condition == "clodronate"),
                     n_sham = sum(rsi_valid$condition == "sham"), sigma_reference = sigma_ref)

md <- md[match(colnames(sce), cell_id)]
stopifnot(identical(md$cell_id, colnames(sce)))
colData(sce) <- S4Vectors::DataFrame(as.data.frame(md), row.names = md$cell_id)
metadata(sce)$annotation_markers <- markers
metadata(sce)$state_programs <- programs
metadata(sce)$ambient_RNA_status <- "not estimable: only Cell Ranger filtered matrices provided"
metadata(sce)$qc_thresholds <- rbindlist(thresholds)
saveRDS(sce, object_out, compress = FALSE)

fwrite(rbindlist(thresholds), file.path(outdir, "qc_thresholds.tsv"), sep = "\t")
fwrite(rbindlist(qc_rows), file.path(outdir, "qc_summary.tsv"), sep = "\t")
fwrite(md, file.path(outdir, "cell_metadata.tsv.gz"), sep = "\t")
fwrite(sample_phase, file.path(outdir, "sample_celltype_phase.tsv"), sep = "\t")
fwrite(rsi, file.path(outdir, "sample_RSI.tsv"), sep = "\t")
fwrite(effect, file.path(outdir, "failure_effect.tsv"), sep = "\t")

p <- ggplot(rsi, aes(condition, RSI, colour = condition)) +
  geom_point(size = 2.3, position = position_jitter(width = 0.05)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, colour = "black") +
  scale_colour_manual(values = c(sham = "#2878B5", clodronate = "#D1495B")) +
  labs(x = NULL, y = "Regenerative Synchrony Index",
       subtitle = sprintf("Exact two-sided permutation P = %.3f", perm_p_two)) +
  theme_classic(base_size = 10) + theme(legend.position = "none")
ggsave(file.path(outdir, "macrophage_depletion_RSI.pdf"), p, width = 3.8, height = 3.5)

cat("GSE291782 complete:", ncol(sce), "QC singlets; RSI effect =",
    signif(observed_effect, 4), "P =", signif(perm_p_two, 4), "\n")

  },
  "hypothyroid" = function() {
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE311479")
outdir <- "results/GSE311479_failure"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(input_dir, pattern = "^[^._].*Normalized_counts.*txt\\.gz$", full.names = TRUE)
stopifnot(length(files) == 6L)

markers <- list(
  MuSC_myogenic = c("Pax7", "Calcr", "Sdc4", "Vcam1", "Myod1", "Myog", "Des"),
  myonucleus = c("Ttn", "Myh1", "Myh2", "Myh4", "Ckm", "Acta1"),
  FAP = c("Pdgfra", "Col1a1", "Col3a1", "Dcn", "Lum"),
  macrophage = c("Ptprc", "Lyz2", "Adgre1", "Cd68", "C1qa", "C1qb"),
  neutrophil = c("S100a8", "S100a9", "Retnlg", "Camp", "Ly6g"),
  T_NK = c("Cd3d", "Cd3e", "Trbc1", "Nkg7", "Klrd1"),
  endothelial = c("Pecam1", "Kdr", "Emcn", "Cdh5", "Klf2"),
  pericyte = c("Rgs5", "Cspg4", "Pdgfrb", "Notch3", "Des")
)
programs <- list(
  Activation = c("Fos", "Jun", "Junb", "Egr1", "Myod1"),
  Proliferation = c("Mki67", "Top2a", "Ccnb1", "Pcna", "Ube2c"),
  Commitment = c("Myod1", "Myf5", "Des", "Cdkn1c"),
  Fusion = c("Myog", "Mymk", "Mymx", "Tnnt1", "Ckm"),
  Remodeling = c("Col1a1", "Col3a1", "Fn1", "Mmp14", "Tgfb1"),
  Resolution = c("Mrc1", "C1qa", "C1qb", "Il10", "Apoe")
)
needed <- unique(c(unlist(markers), unlist(programs)))

read_selected <- function(f) {
  z <- fread(f, showProgress = FALSE)
  gene <- z[[1]]
  cell_cols <- names(z)[14:ncol(z)]
  hit <- which(gene %in% needed)
  m <- as.matrix(z[hit, ..cell_cols])
  storage.mode(m) <- "double"
  rownames(m) <- gene[hit]
  sample_id <- sub("\\.txt\\.gz$", "", basename(f))
  colnames(m) <- paste(sample_id, cell_cols, sep = "::")
  rm(z); gc()
  m
}

mats <- lapply(files, read_selected)
all_genes <- unique(unlist(lapply(mats, rownames)))
aligned <- lapply(mats, function(m) {
  out <- matrix(0, nrow = length(all_genes), ncol = ncol(m),
                dimnames = list(all_genes, colnames(m)))
  out[rownames(m), ] <- m
  out
})
expr <- do.call(cbind, aligned)
rm(mats, aligned); gc()

cell_id <- colnames(expr)
sample_id <- sub("::.*", "", cell_id)
condition <- ifelse(grepl("CTRL", sample_id), "control", "hypothyroid")
time_days <- ifelse(grepl("_NI$", sample_id), 0,
                    ifelse(grepl("_05d$", sample_id), 5, 14))

# The files contain author-normalized values, not raw counts. They support
# expression-state scoring but not count-based QC, ambient correction or doublet calling.
feature_scaled <- scale(t(expr))
feature_scaled[!is.finite(feature_scaled)] <- 0
pca <- prcomp(feature_scaled, center = FALSE, scale. = FALSE, rank. = 12)
set.seed(20260824)
cluster <- kmeans(pca$x[, seq_len(min(12, ncol(pca$x))), drop = FALSE],
                  centers = 18, nstart = 20, iter.max = 100)$cluster

score_module <- function(gs) {
  hit <- intersect(gs, rownames(expr))
  if (!length(hit)) return(rep(NA_real_, ncol(expr)))
  colMeans(expr[hit, , drop = FALSE])
}
md <- data.table(cell_id, sample_id, condition, time_days, cluster_id = as.character(cluster))
for (nm in names(markers)) md[[paste0("marker_", nm)]] <- score_module(markers[[nm]])
for (nm in names(programs)) md[[paste0("state_score_", nm)]] <- score_module(programs[[nm]])

marker_cols <- paste0("marker_", names(markers))
cl <- md[, lapply(.SD, mean, na.rm = TRUE), by = cluster_id, .SDcols = marker_cols]
M <- as.matrix(cl[, ..marker_cols])
top <- max.col(M, ties.method = "first")
ordered <- t(apply(M, 1, sort, decreasing = TRUE))
cl[, coarse_cell_type := names(markers)[top]]
cl[, annotation_margin := ordered[, 1] - ordered[, 2]]
cl[!is.finite(annotation_margin) | annotation_margin < 0.08,
   coarse_cell_type := "Unknown_or_ambiguous"]
md <- merge(md, cl[, .(cluster_id, coarse_cell_type, annotation_margin)],
            by = "cluster_id", sort = FALSE)

program_cols <- paste0("state_score_", names(programs))
sample_phase <- md[coarse_cell_type != "Unknown_or_ambiguous",
  c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
  by = .(sample_id, condition, time_days, coarse_cell_type), .SDcols = program_cols]
for (pc in program_cols) {
  mu <- mean(sample_phase[[pc]], na.rm = TRUE); ss <- sd(sample_phase[[pc]], na.rm = TRUE)
  sample_phase[[pc]] <- (sample_phase[[pc]] - mu) / ifelse(ss > 0, ss, 1)
}
centers <- c(0.12, 0.28, 0.45, 0.64, 0.80, 0.94)
z <- as.matrix(sample_phase[, ..program_cols])
z <- sweep(z, 1, apply(z, 1, max), "-") / 0.8
pp <- exp(z); pp <- pp / rowSums(pp)
sample_phase[, phase := as.numeric(pp %*% centers)]
# Safer explicit same-time, same-lineage control lookup.
ctrl <- sample_phase[condition == "control", .(time_days, coarse_cell_type, expected = phase)]
sample_phase <- merge(sample_phase, ctrl,
                      by = c("time_days", "coarse_cell_type"), all.x = TRUE)
sample_phase[, deviation := phase - expected]
sample_phase[, reliability := pmin(1, sqrt(n_cells / 100))]

hypo_rsi <- sample_phase[condition == "hypothyroid" & n_cells >= 20 & is.finite(deviation),
  .(RSI_descriptive = regenerative_synchrony_index(deviation, reliability,
                                                    sigma_reference = 0.10),
    dispersion = weighted_var(deviation, reliability), effective_cell_types = .N),
  by = .(sample_id, time_days)]

program_contrasts <- merge(
  sample_phase[condition == "hypothyroid"],
  sample_phase[condition == "control"],
  by = c("time_days", "coarse_cell_type"), suffixes = c("_hypo", "_control"))
for (pc in program_cols) {
  program_contrasts[[paste0(pc, "_difference")]] <-
    program_contrasts[[paste0(pc, "_hypo")]] - program_contrasts[[paste0(pc, "_control")]]
}

fwrite(md, file.path(outdir, "cell_metadata_selected_programs.tsv.gz"), sep = "\t")
fwrite(sample_phase, file.path(outdir, "sample_celltype_phase.tsv"), sep = "\t")
fwrite(hypo_rsi, file.path(outdir, "hypothyroid_RSI_descriptive.tsv"), sep = "\t")
fwrite(program_contrasts, file.path(outdir, "program_strength_contrasts.tsv"), sep = "\t")
writeLines(c(
  "inference_status: descriptive_only",
  "reason: each condition-time expression matrix pools three mice into one library",
  "input_status: author-normalized values; raw-count QC/scDblFinder/ambient correction unavailable",
  "RSI_status: same-time displacement from control, not a population-level estimate"
), file.path(outdir, "analysis_manifest.txt"))

p <- ggplot(sample_phase[n_cells >= 20],
            aes(coarse_cell_type, deviation, fill = condition)) +
  geom_col(position = "dodge") + facet_wrap(~time_days, labeller = label_both) +
  scale_fill_manual(values = c(control = "#2878B5", hypothyroid = "#D1495B")) +
  labs(x = NULL, y = "Phase displacement from same-time control", fill = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(outdir, "hypothyroid_phase_displacement.pdf"), p,
       width = 7.2, height = 4)

cat("GSE311479 descriptive failure analysis complete:", ncol(expr), "cells.\n")

  },
  "notexin_extension" = function() {
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(BiocParallel)
  library(data.table)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE159500")
outdir <- "results/GSE159500_extension"
object_out <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "processed/GSE159500_qc_singlet_sce.rds")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_out), recursive = TRUE, showWarnings = FALSE)

prefixes <- sub("_soupx_matrix.mtx.gz$", "",
                list.files(input_dir, pattern = "^[^._].*_soupx_matrix.mtx.gz$",
                           full.names = TRUE))
stopifnot(length(prefixes) == 2L)

map_author_type <- function(x) {
  fcase(
    grepl("Satellite|MuSC|Myoblast", x, ignore.case = TRUE), "MuSC_myogenic",
    grepl("Myonuc", x, ignore.case = TRUE), "myonucleus",
    grepl("FAP|Fibro|Tenocyte", x, ignore.case = TRUE), "FAP",
    grepl("Macro|Mono", x, ignore.case = TRUE), "macrophage",
    grepl("Neutro", x, ignore.case = TRUE), "neutrophil",
    grepl("T Cell|NK", x, ignore.case = TRUE), "T_NK",
    grepl("Endothel", x, ignore.case = TRUE), "endothelial",
    grepl("Pericyte|Smooth Muscle", x, ignore.case = TRUE), "pericyte",
    default = "other_major"
  )
}

set.seed(20260824)
objects <- list(); summaries <- list()
for (prefix in prefixes) {
  mat <- as(readMM(gzfile(paste0(prefix, "_soupx_matrix.mtx.gz"))), "dgCMatrix")
  genes <- fread(paste0(prefix, "_soupx_genes.tsv.gz"), header = FALSE)[[1]]
  bc <- fread(paste0(prefix, "_soupx_barcodes.tsv.gz"), header = FALSE)[[1]]
  rownames(mat) <- make.unique(genes); colnames(mat) <- bc
  meta_file <- sub("_soupx$", "", prefix)
  meta_file <- paste0(meta_file, "_meta.csv.gz")
  meta <- fread(meta_file)
  meta_barcode <- sub("^_", "", meta[[1]])
  idx <- match(sub("^_", "", bc), meta_barcode)
  # SoupX matrices include additional corrected droplets; the author metadata
  # defines the retained biological cells and is therefore the inclusion set.
  keep <- !is.na(idx)
  if (sum(keep) < 0.95 * nrow(meta)) stop("Fewer than 95% of annotated cells matched")
  mat <- mat[, keep, drop = FALSE]
  bc <- bc[keep]
  idx <- idx[keep]
  colnames(mat) <- bc
  sample_id <- unique(meta$sample)[1]
  sce <- SingleCellExperiment(list(counts = mat))
  sce$sample_id <- sample_id
  sce$time_days <- unique(meta$injury.days)[1]
  sce$author_cell_type <- meta$cell.types[idx]
  sce$coarse_cell_type <- map_author_type(sce$author_cell_type)
  mt <- grepl("^mt-", rownames(sce), ignore.case = TRUE)
  qc <- perCellQCMetrics(sce, subsets = list(Mt = mt))
  low_detected <- max(200, median(qc$detected) - 3 * mad(qc$detected))
  low_sum <- max(500, median(qc$sum) - 3 * mad(qc$sum))
  high_detected <- median(qc$detected) + 5 * mad(qc$detected)
  high_sum <- median(qc$sum) + 5 * mad(qc$sum)
  high_mt <- min(40, median(qc$subsets_Mt_percent) + 3 * mad(qc$subsets_Mt_percent))
  pass <- qc$detected >= low_detected & qc$sum >= low_sum &
    qc$detected <= high_detected & qc$sum <= high_sum & qc$subsets_Mt_percent <= high_mt
  colData(sce) <- cbind(colData(sce), qc)
  sce$passes_core_qc <- pass
  cls <- rep(NA_character_, ncol(sce)); score <- rep(NA_real_, ncol(sce))
  called <- scDblFinder(sce[, pass], samples = rep(sample_id, sum(pass)),
                        BPPARAM = SerialParam(), verbose = FALSE)
  cls[pass] <- as.character(called$scDblFinder.class)
  score[pass] <- called$scDblFinder.score
  sce$scDblFinder.class <- cls; sce$scDblFinder.score <- score
  sce$passes_qc_and_singlet <- pass & cls == "singlet"
  summaries[[sample_id]] <- data.table(sample_id, input_cells = ncol(sce),
    core_qc_pass = sum(pass), doublets = sum(cls == "doublet", na.rm = TRUE),
    final_singlets = sum(sce$passes_qc_and_singlet, na.rm = TRUE),
    low_detected, low_sum, high_detected, high_sum, high_mt)
  objects[[sample_id]] <- sce[, sce$passes_qc_and_singlet %in% TRUE]
}

sce <- do.call(cbind, objects)
sizeFactors(sce) <- NULL
sce <- logNormCounts(sce)
programs <- list(
  Activation = c("Fos", "Jun", "Junb", "Egr1", "Myod1"),
  Proliferation = c("Mki67", "Top2a", "Ccnb1", "Pcna", "Ube2c"),
  Commitment = c("Myod1", "Myf5", "Des", "Cdkn1c"),
  Fusion = c("Myog", "Mymk", "Mymx", "Tnnt1", "Ckm"),
  Remodeling = c("Col1a1", "Col3a1", "Fn1", "Mmp14", "Tgfb1"),
  Resolution = c("Mrc1", "C1qa", "C1qb", "Il10", "Apoe")
)
for (nm in names(programs)) {
  hit <- intersect(programs[[nm]], rownames(sce))
  colData(sce)[[paste0("state_score_", nm)]] <-
    if (length(hit)) Matrix::colMeans(logcounts(sce)[hit, , drop = FALSE]) else NA_real_
}

md <- as.data.table(as.data.frame(colData(sce)), keep.rownames = "cell_id")
program_cols <- paste0("state_score_", names(programs))
sample_state <- md[, c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
                   by = .(sample_id, time_days, coarse_cell_type), .SDcols = program_cols]
fwrite(rbindlist(summaries), file.path(outdir, "qc_summary.tsv"), sep = "\t")
fwrite(sample_state, file.path(outdir, "sample_celltype_states.tsv"), sep = "\t")
fwrite(md, file.path(outdir, "cell_metadata.tsv.gz"), sep = "\t")
metadata(sce)$ambient_correction <- "author-provided SoupX matrix"
saveRDS(sce, object_out, compress = FALSE)
cat("GSE159500 complete:", ncol(sce), "QC singlets.\n")

  },
  "barium_chloride" = function() {
suppressPackageStartupMessages({
  library(DropletUtils)
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(scater)
  library(bluster)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE205690")
outdir <- "results/GSE205690_bacl2"
object_out <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "processed/GSE205690_qc_annotated_h5sce")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_out), recursive = TRUE, showWarnings = FALSE)
files <- list.files(input_dir, pattern = "^[^._].*raw_feature_bc_matrix.h5$", full.names = TRUE)
stopifnot(length(files) == 24L)

set.seed(20260824)
objects <- list(); qc_list <- list(); threshold_list <- list()
file_unit <- vapply(files, function(f) strsplit(sub("_raw_feature_bc_matrix.h5$", "",
                                                   basename(f)), "_")[[1]][2], character(1))
file_groups <- split(files, file_unit)
stopifnot(length(file_groups) == 6L, all(lengths(file_groups) == 4L))
for (bio_unit in names(file_groups)) {
  group_files <- file_groups[[bio_unit]]
  raw_parts <- lapply(group_files, function(f) read10xCounts(f, type = "HDF5", col.names = TRUE))
  stopifnot(all(vapply(raw_parts[-1], function(z)
    identical(colnames(z), colnames(raw_parts[[1]])), logical(1))))
  raw_counts <- Reduce(`+`, lapply(raw_parts, function(z) as(counts(z), "dgCMatrix")))
  raw <- SingleCellExperiment(list(counts = raw_counts),
                              rowData = rowData(raw_parts[[1]]))
  rownames(raw) <- make.unique(rowData(raw)$Symbol)
  totals <- Matrix::colSums(counts(raw))
  br <- barcodeRanks(counts(raw), lower = 100)
  inflection <- metadata(br)$inflection
  cell_threshold <- max(100, ifelse(is.finite(inflection), inflection, 100))
  called <- which(totals >= cell_threshold)
  if (length(called) < 50L) stop("Fewer than 50 cell-containing droplets in animal ", bio_unit)

  empty <- which(totals > 0 & totals <= 50)
  if (length(empty) > 50000L) empty <- sample(empty, 50000L)
  ambient <- Matrix::rowSums(counts(raw)[, empty, drop = FALSE])
  top_ambient <- order(ambient, decreasing = TRUE)[seq_len(min(20L, sum(ambient > 0)))]

  sce <- raw[, called]
  sample_tag <- sub("_[1-4]_raw_feature_bc_matrix.h5$", "", basename(group_files[[1]]))
  time_days <- if (grepl("YoungUI", sample_tag)) 0 else if (grepl("Day4", sample_tag)) 4 else 7
  capture_id <- paste(sub("_raw_feature_bc_matrix.h5$", "", basename(group_files)), collapse = ";")
  colnames(sce) <- paste0("animal_", bio_unit, "::", colnames(sce))
  sce$sample_id <- paste0("animal_", bio_unit)
  sce$capture_id <- capture_id
  sce$time_days <- time_days
  qc <- perCellQCMetrics(sce, subsets = list(Mt = grepl("^mt-", rownames(sce), ignore.case = TRUE)))
  qc$ambient_top20_percent <- 100 * Matrix::colSums(counts(sce)[top_ambient, , drop = FALSE]) /
    pmax(qc$sum, 1)
  low_detected <- max(150, median(qc$detected) - 3 * mad(qc$detected))
  low_sum <- max(300, median(qc$sum) - 3 * mad(qc$sum))
  high_detected <- median(qc$detected) + 5 * mad(qc$detected)
  high_sum <- median(qc$sum) + 5 * mad(qc$sum)
  high_mt <- min(40, median(qc$subsets_Mt_percent) + 3 * mad(qc$subsets_Mt_percent))
  pass <- qc$detected >= low_detected & qc$sum >= low_sum &
    qc$detected <= high_detected & qc$sum <= high_sum & qc$subsets_Mt_percent <= high_mt
  colData(sce) <- cbind(colData(sce), qc)
  sce$passes_core_qc <- pass
  dbl <- scDblFinder(sce[, pass], samples = rep(paste0("animal_", bio_unit), sum(pass)),
                     BPPARAM = SerialParam(), verbose = FALSE)
  cls <- rep(NA_character_, ncol(sce)); score <- rep(NA_real_, ncol(sce))
  cls[pass] <- as.character(dbl$scDblFinder.class); score[pass] <- dbl$scDblFinder.score
  sce$scDblFinder.class <- cls; sce$scDblFinder.score <- score
  sce$passes_qc_and_singlet <- pass & cls == "singlet"
  objects[[bio_unit]] <- sce[, sce$passes_qc_and_singlet %in% TRUE]
  qc_list[[bio_unit]] <- data.table(capture_id, sample_id = unique(sce$sample_id),
    time_days, technical_files_summed = length(group_files), raw_droplets = ncol(raw),
    called_droplets = length(called),
    core_qc_pass = sum(pass), doublets = sum(cls == "doublet", na.rm = TRUE),
    final_singlets = sum(sce$passes_qc_and_singlet, na.rm = TRUE),
    median_ambient_top20_percent = median(qc$ambient_top20_percent, na.rm = TRUE))
  threshold_list[[bio_unit]] <- data.table(capture_id, sample_id = unique(sce$sample_id),
                                              cell_threshold, low_detected,
                                              low_sum, high_detected, high_sum, high_mt)
  rm(raw_parts, raw_counts, raw, sce, dbl); gc()
}

sce <- do.call(cbind, objects)
rm(objects); gc()
sizeFactors(sce) <- NULL
sce <- logNormCounts(sce)
gene_var <- DelayedMatrixStats::rowVars(logcounts(sce), na.rm = TRUE)
hvg <- head(order(gene_var, decreasing = TRUE, na.last = NA), 3000L)
sce <- runPCA(sce, subset_row = hvg, ncomponents = 25,
              BSPARAM = BiocSingular::IrlbaParam())
sce$cluster_id <- as.character(clusterRows(reducedDim(sce, "PCA"),
  NNGraphParam(k = 15, cluster.fun = "louvain")))

markers <- list(
  MuSC_myogenic = c("Pax7", "Calcr", "Sdc4", "Vcam1", "Myod1", "Myog", "Des"),
  myonucleus = c("Ttn", "Myh1", "Myh2", "Myh4", "Ckm", "Acta1"),
  FAP = c("Pdgfra", "Col1a1", "Col3a1", "Dcn", "Lum"),
  macrophage = c("Ptprc", "Lyz2", "Adgre1", "Cd68", "C1qa", "C1qb"),
  neutrophil = c("S100a8", "S100a9", "Retnlg", "Camp", "Ly6g"),
  T_NK = c("Cd3d", "Cd3e", "Trbc1", "Nkg7", "Klrd1"),
  endothelial = c("Pecam1", "Kdr", "Emcn", "Cdh5", "Klf2"),
  pericyte = c("Rgs5", "Cspg4", "Pdgfrb", "Notch3", "Des")
)
programs <- list(
  Activation = c("Fos", "Jun", "Junb", "Egr1", "Myod1"),
  Proliferation = c("Mki67", "Top2a", "Ccnb1", "Pcna", "Ube2c"),
  Commitment = c("Myod1", "Myf5", "Des", "Cdkn1c"),
  Fusion = c("Myog", "Mymk", "Mymx", "Tnnt1", "Ckm"),
  Remodeling = c("Col1a1", "Col3a1", "Fn1", "Mmp14", "Tgfb1"),
  Resolution = c("Mrc1", "C1qa", "C1qb", "Il10", "Apoe")
)
score_module <- function(gs) {
  hit <- intersect(gs, rownames(sce))
  if (!length(hit)) return(rep(NA_real_, ncol(sce)))
  Matrix::colMeans(logcounts(sce)[hit, , drop = FALSE])
}
for (nm in names(markers)) colData(sce)[[paste0("marker_", nm)]] <- score_module(markers[[nm]])
for (nm in names(programs)) colData(sce)[[paste0("state_score_", nm)]] <- score_module(programs[[nm]])
md0 <- as.data.frame(colData(sce))
md <- as.data.table(md0, keep.rownames = "cell_id")
marker_cols <- paste0("marker_", names(markers))
cl <- md[, lapply(.SD, mean, na.rm = TRUE), by = cluster_id, .SDcols = marker_cols]
M <- as.matrix(cl[, ..marker_cols]); top <- max.col(M, ties.method = "first")
ordered <- t(apply(M, 1, sort, decreasing = TRUE))
cl[, coarse_cell_type := names(markers)[top]]
cl[, annotation_margin := ordered[, 1] - ordered[, 2]]
cl[!is.finite(annotation_margin) | annotation_margin < 0.08,
   coarse_cell_type := "Unknown_or_ambiguous"]
md <- merge(md, cl[, .(cluster_id, coarse_cell_type, annotation_margin)],
            by = "cluster_id", sort = FALSE)
md <- md[match(colnames(sce), cell_id)]
stopifnot(identical(md$cell_id, colnames(sce)))
colData(sce) <- S4Vectors::DataFrame(as.data.frame(md), row.names = md$cell_id)

program_cols <- paste0("state_score_", names(programs))
sample_state <- md[coarse_cell_type != "Unknown_or_ambiguous",
  c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
  by = .(sample_id, time_days, coarse_cell_type), .SDcols = program_cols]

fwrite(rbindlist(qc_list), file.path(outdir, "qc_summary_by_biological_sample.tsv"), sep = "\t")
fwrite(rbindlist(threshold_list), file.path(outdir, "qc_thresholds_by_biological_sample.tsv"), sep = "\t")
fwrite(md, file.path(outdir, "cell_metadata.tsv.gz"), sep = "\t")
fwrite(sample_state, file.path(outdir, "sample_celltype_states.tsv"), sep = "\t")
metadata(sce)$ambient_metric <- "fraction of counts in top 20 genes from <=50-UMI droplets; review signal, not hard filter"
metadata(sce)$technical_recombination <- "four technical count matrices summed by matching barcode within animal before cell calling"
if (!dir.exists(object_out)) {
  HDF5Array::saveHDF5SummarizedExperiment(sce, dir = object_out, replace = FALSE)
} else {
  message("Preserving existing HDF5 checkpoint: ", object_out)
}

p <- ggplot(sample_state, aes(time_days, n_cells, colour = coarse_cell_type,
                              group = interaction(sample_id, coarse_cell_type))) +
  geom_point(position = position_jitter(width = 0.08), size = 1.6) + scale_y_log10() +
  labs(x = "Days after BaCl2 injury", y = "Cells per biological sample", colour = NULL) +
  theme_classic(base_size = 9)
ggsave(file.path(outdir, "celltype_dynamics.pdf"), p, width = 6.5, height = 4)
cat("GSE205690 complete:", ncol(sce), "QC singlets across 6 biological samples.\n")

  },
  "old_notexin" = function() {
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE162172")
outdir <- "results/GSE162172_old"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
zip_files <- list.files(input_dir, pattern = "^[^._].*_soupx.zip$", full.names = TRUE)
stopifnot(length(zip_files) == 21L)

map_author_type <- function(x) {
  fcase(
    grepl("Satellite|MuSC|Myoblast|Myogenic", x, ignore.case = TRUE), "MuSC_myogenic",
    grepl("Myonuc", x, ignore.case = TRUE), "myonucleus",
    grepl("FAP|Fibro|Tenocyte", x, ignore.case = TRUE), "FAP",
    grepl("Macro|Mono", x, ignore.case = TRUE), "macrophage",
    grepl("Neutro", x, ignore.case = TRUE), "neutrophil",
    grepl("T Cell|NK", x, ignore.case = TRUE), "T_NK",
    grepl("Endothel", x, ignore.case = TRUE), "endothelial",
    grepl("Pericyte|Smooth Muscle", x, ignore.case = TRUE), "pericyte",
    default = "other_major"
  )
}

programs <- list(
  Activation = c("Fos", "Jun", "Junb", "Egr1", "Myod1"),
  Proliferation = c("Mki67", "Top2a", "Ccnb1", "Pcna", "Ube2c"),
  Commitment = c("Myod1", "Myf5", "Des", "Cdkn1c"),
  Fusion = c("Myog", "Mymk", "Mymx", "Tnnt1", "Ckm"),
  Remodeling = c("Col1a1", "Col3a1", "Fn1", "Mmp14", "Tgfb1"),
  Resolution = c("Mrc1", "C1qa", "C1qb", "Il10", "Apoe")
)

set.seed(20260824)
qc_all <- list(); meta_all <- list(); state_all <- list()
for (zipf in zip_files) {
  members <- unzip(zipf, list = TRUE)$Name
  mm_name <- members[grepl("/matrix.mtx$", members)]
  gene_name <- members[grepl("/genes.tsv$", members)]
  bc_name <- members[grepl("/barcodes.tsv$", members)]
  stopifnot(length(mm_name) == 1L, length(gene_name) == 1L, length(bc_name) == 1L)
  mat <- as(readMM(unz(zipf, mm_name)), "dgCMatrix")
  genes <- read.delim(unz(zipf, gene_name), header = FALSE, stringsAsFactors = FALSE)[[1]]
  bc <- read.delim(unz(zipf, bc_name), header = FALSE, stringsAsFactors = FALSE)[[1]]
  rownames(mat) <- make.unique(genes); colnames(mat) <- bc

  sample_stub <- sub("_soupx.zip$", "", basename(zipf))
  meta_file <- file.path(input_dir, paste0(sample_stub, "_meta.csv.gz"))
  meta <- fread(meta_file)
  sample_id <- unique(meta$sample)[1]
  meta_bc <- sub(paste0("^", sample_id, "_"), "", meta[[1]])
  idx <- match(bc, meta_bc)
  keep <- !is.na(idx)
  if (sum(keep) < 0.95 * nrow(meta)) stop("Annotated-cell match below 95% for ", sample_id)
  mat <- mat[, keep, drop = FALSE]; bc <- bc[keep]; idx <- idx[keep]
  colnames(mat) <- paste(sample_id, bc, sep = "::")

  sce <- SingleCellExperiment(list(counts = mat))
  sce$sample_id <- sample_id
  sce$time_days <- unique(meta$injury.days)[1]
  sce$replicate <- sub(".*_", "", sample_id)
  sce$author_doublet_call <- as.character(meta$DF.individual[idx])
  sce$author_noisy <- as.logical(meta$is.noisy[idx])
  sce$author_cell_type <- as.character(meta$deconIDs[idx])
  sce$coarse_cell_type <- map_author_type(sce$author_cell_type)

  qc <- perCellQCMetrics(sce, subsets = list(Mt = grepl("^mt-", rownames(sce), ignore.case = TRUE)))
  low_detected <- max(150, median(qc$detected) - 3 * mad(qc$detected))
  low_sum <- max(300, median(qc$sum) - 3 * mad(qc$sum))
  high_detected <- median(qc$detected) + 5 * mad(qc$detected)
  high_sum <- median(qc$sum) + 5 * mad(qc$sum)
  high_mt <- min(40, median(qc$subsets_Mt_percent) + 3 * mad(qc$subsets_Mt_percent))
  pass <- qc$detected >= low_detected & qc$sum >= low_sum &
    qc$detected <= high_detected & qc$sum <= high_sum & qc$subsets_Mt_percent <= high_mt
  colData(sce) <- cbind(colData(sce), qc)
  sce$passes_core_qc <- pass
  dbl <- scDblFinder(sce[, pass], samples = rep(sample_id, sum(pass)),
                     BPPARAM = SerialParam(), verbose = FALSE)
  cls <- rep(NA_character_, ncol(sce)); score <- rep(NA_real_, ncol(sce))
  cls[pass] <- as.character(dbl$scDblFinder.class); score[pass] <- dbl$scDblFinder.score
  sce$scDblFinder.class <- cls; sce$scDblFinder.score <- score
  sce$passes_qc_and_singlet <- pass & cls == "singlet"
  final <- sce[, sce$passes_qc_and_singlet %in% TRUE]
  sizeFactors(final) <- NULL
  final <- logNormCounts(final)
  for (nm in names(programs)) {
    hit <- intersect(programs[[nm]], rownames(final))
    colData(final)[[paste0("state_score_", nm)]] <-
      if (length(hit)) Matrix::colMeans(logcounts(final)[hit, , drop = FALSE]) else NA_real_
  }
  md <- as.data.table(as.data.frame(colData(final)), keep.rownames = "cell_id")
  program_cols <- paste0("state_score_", names(programs))
  state <- md[, c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
              by = .(sample_id, time_days, replicate, coarse_cell_type),
              .SDcols = program_cols]
  qc_all[[sample_id]] <- data.table(sample_id, time_days = unique(sce$time_days),
    input_annotated_cells = ncol(sce), core_qc_pass = sum(pass),
    scDblFinder_doublets = sum(cls == "doublet", na.rm = TRUE),
    final_singlets = ncol(final), author_doublets = sum(sce$author_doublet_call != "Singlet", na.rm = TRUE),
    low_detected, low_sum, high_detected, high_sum, high_mt)
  meta_all[[sample_id]] <- md
  state_all[[sample_id]] <- state
  rm(mat, sce, final, dbl, md); gc()
}

qc_tab <- rbindlist(qc_all); cell_meta <- rbindlist(meta_all, fill = TRUE)
sample_state <- rbindlist(state_all, fill = TRUE)
fwrite(qc_tab, file.path(outdir, "qc_summary.tsv"), sep = "\t")
fwrite(cell_meta, file.path(outdir, "cell_inclusion_and_scores.tsv.gz"), sep = "\t")
fwrite(sample_state, file.path(outdir, "sample_celltype_states.tsv"), sep = "\t")
writeLines(c(
  "counts: author SoupX-corrected matrices retained in original zip archives",
  "filtered_object: cell inclusion and state-score manifest; matrices are not duplicated to conserve disk",
  "doublet_detection: scDblFinder rerun per biological sample",
  "inference_unit: one mouse per sample"
), file.path(outdir, "analysis_manifest.txt"))

p <- ggplot(sample_state, aes(time_days, n_cells, colour = coarse_cell_type)) +
  geom_point(position = position_jitter(width = 0.08), alpha = 0.75) +
  stat_summary(fun = median, geom = "line", aes(group = coarse_cell_type)) +
  scale_y_log10() + labs(x = "Days after notexin injury", y = "Cells per mouse", colour = NULL) +
  theme_classic(base_size = 9)
ggsave(file.path(outdir, "old_celltype_dynamics.pdf"), p, width = 6.5, height = 4)
cat("GSE162172 complete:", sum(qc_tab$final_singlets), "QC singlets across",
    nrow(qc_tab), "mice.\n")

  },
  "human_exercise" = function() {
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(scater)
  library(bluster)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE214544")
outdir <- "results/GSE214544_human"
object_out <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "processed/GSE214544_qc_annotated_h5sce")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_out), recursive = TRUE, showWarnings = FALSE)

# GSM6611299/1300 are fully contained in the deeper-sequenced GSM6611301/1302
# barcode sets for donor 3 and are excluded to prevent duplicate cells.
keep_gsm <- c("GSM6611295", "GSM6611296", "GSM6611297", "GSM6611298",
              "GSM6611301", "GSM6611302")
matrix_files <- list.files(input_dir, pattern = "^[^._].*_matrix.mtx.gz$", full.names = TRUE)
matrix_files <- matrix_files[vapply(matrix_files, function(f)
  any(startsWith(basename(f), keep_gsm)), logical(1))]
stopifnot(length(matrix_files) == 6L)

sample_design <- data.table(
  gsm = keep_gsm,
  donor = rep(c("donor_1", "donor_2", "donor_3"), each = 2),
  condition = rep(c("pre", "post_3h"), 3)
)

read_bundle <- function(f) {
  prefix <- sub("_matrix.mtx.gz$", "", f)
  feat <- fread(paste0(prefix, "_features.tsv.gz"), header = FALSE)
  bc <- fread(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)[[1]]
  m <- as(readMM(gzfile(f)), "dgCMatrix")
  rownames(m) <- make.unique(feat[[2]])
  gsm_id <- sub("_.*", "", basename(prefix))
  colnames(m) <- paste(gsm_id, bc, sep = "::")
  design <- sample_design[gsm == gsm_id]
  sce <- SingleCellExperiment(list(counts = m))
  sce$sample_id <- gsm_id; sce$donor <- design$donor; sce$condition <- design$condition
  sce
}

set.seed(20260824)
objects <- lapply(matrix_files, read_bundle)
qc_rows <- list(); threshold_rows <- list()
for (i in seq_along(objects)) {
  sce <- objects[[i]]; sid <- unique(sce$sample_id)
  qc <- perCellQCMetrics(sce, subsets = list(Mt = grepl("^MT-", rownames(sce))))
  low_detected <- max(200, median(qc$detected) - 3 * mad(qc$detected))
  low_sum <- max(500, median(qc$sum) - 3 * mad(qc$sum))
  high_detected <- median(qc$detected) + 5 * mad(qc$detected)
  high_sum <- median(qc$sum) + 5 * mad(qc$sum)
  high_mt <- min(40, median(qc$subsets_Mt_percent) + 3 * mad(qc$subsets_Mt_percent))
  pass <- qc$detected >= low_detected & qc$sum >= low_sum &
    qc$detected <= high_detected & qc$sum <= high_sum & qc$subsets_Mt_percent <= high_mt
  colData(sce) <- cbind(colData(sce), qc); sce$passes_core_qc <- pass
  dbl <- scDblFinder(sce[, pass], samples = rep(sid, sum(pass)),
                     BPPARAM = SerialParam(), verbose = FALSE)
  cls <- rep(NA_character_, ncol(sce)); score <- rep(NA_real_, ncol(sce))
  cls[pass] <- as.character(dbl$scDblFinder.class); score[pass] <- dbl$scDblFinder.score
  sce$scDblFinder.class <- cls; sce$scDblFinder.score <- score
  sce$passes_qc_and_singlet <- pass & cls == "singlet"
  qc_rows[[sid]] <- data.table(sample_id = sid, donor = unique(sce$donor),
    condition = unique(sce$condition), input_cells = ncol(sce), core_qc_pass = sum(pass),
    doublets = sum(cls == "doublet", na.rm = TRUE),
    final_singlets = sum(sce$passes_qc_and_singlet, na.rm = TRUE))
  threshold_rows[[sid]] <- data.table(sample_id = sid, low_detected, low_sum,
                                       high_detected, high_sum, high_mt)
  objects[[i]] <- sce[, sce$passes_qc_and_singlet %in% TRUE]
}

sce <- do.call(cbind, objects); rm(objects); gc()
sizeFactors(sce) <- NULL; sce <- logNormCounts(sce)
gene_var <- DelayedMatrixStats::rowVars(logcounts(sce), na.rm = TRUE)
hvg <- head(order(gene_var, decreasing = TRUE, na.last = NA), 3000L)
sce <- runPCA(sce, subset_row = hvg, ncomponents = 30,
              BSPARAM = BiocSingular::IrlbaParam())
sce$cluster_id <- as.character(clusterRows(reducedDim(sce, "PCA"),
  NNGraphParam(k = 20, cluster.fun = "louvain")))

markers <- list(
  MuSC_myogenic = c("PAX7", "CALCR", "SDC4", "VCAM1", "MYOD1", "MYOG", "DES"),
  myonucleus = c("TTN", "MYH1", "MYH2", "MYH4", "CKM", "ACTA1"),
  FAP = c("PDGFRA", "COL1A1", "COL3A1", "DCN", "LUM"),
  macrophage = c("PTPRC", "LYZ", "ADGRE1", "CD68", "C1QA", "C1QB"),
  neutrophil = c("S100A8", "S100A9", "FCGR3B", "CSF3R"),
  T_NK = c("CD3D", "CD3E", "TRBC1", "NKG7", "KLRD1"),
  endothelial = c("PECAM1", "KDR", "EMCN", "CDH5", "KLF2"),
  pericyte = c("RGS5", "CSPG4", "PDGFRB", "NOTCH3", "DES")
)
programs <- list(
  Activation = c("FOS", "JUN", "JUNB", "EGR1", "MYOD1"),
  Proliferation = c("MKI67", "TOP2A", "CCNB1", "PCNA", "UBE2C"),
  Commitment = c("MYOD1", "MYF5", "DES", "CDKN1C"),
  Fusion = c("MYOG", "MYMK", "MYMX", "TNNT1", "CKM"),
  Remodeling = c("COL1A1", "COL3A1", "FN1", "MMP14", "TGFB1"),
  Resolution = c("MRC1", "C1QA", "C1QB", "IL10", "APOE")
)
score_module <- function(gs) {
  hit <- intersect(gs, rownames(sce))
  if (!length(hit)) return(rep(NA_real_, ncol(sce)))
  Matrix::colMeans(logcounts(sce)[hit, , drop = FALSE])
}
for (nm in names(markers)) colData(sce)[[paste0("marker_", nm)]] <- score_module(markers[[nm]])
for (nm in names(programs)) colData(sce)[[paste0("state_score_", nm)]] <- score_module(programs[[nm]])
md <- as.data.table(as.data.frame(colData(sce)), keep.rownames = "cell_id")
marker_cols <- paste0("marker_", names(markers))
cl <- md[, lapply(.SD, mean, na.rm = TRUE), by = cluster_id, .SDcols = marker_cols]
M <- as.matrix(cl[, ..marker_cols]); top <- max.col(M, ties.method = "first")
ordered <- t(apply(M, 1, sort, decreasing = TRUE))
cl[, coarse_cell_type := names(markers)[top]]
cl[, annotation_margin := ordered[, 1] - ordered[, 2]]
cl[!is.finite(annotation_margin) | annotation_margin < 0.08,
   coarse_cell_type := "Unknown_or_ambiguous"]
md <- merge(md, cl[, .(cluster_id, coarse_cell_type, annotation_margin)],
            by = "cluster_id", sort = FALSE)
md <- md[match(colnames(sce), cell_id)]
colData(sce) <- S4Vectors::DataFrame(as.data.frame(md), row.names = md$cell_id)

program_cols <- paste0("state_score_", names(programs))
sample_phase <- md[coarse_cell_type != "Unknown_or_ambiguous",
  c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
  by = .(sample_id, donor, condition, coarse_cell_type), .SDcols = program_cols]
for (pc in program_cols) {
  mu <- mean(sample_phase[[pc]], na.rm = TRUE); ss <- sd(sample_phase[[pc]], na.rm = TRUE)
  sample_phase[[pc]] <- (sample_phase[[pc]] - mu) / ifelse(ss > 0, ss, 1)
}
centers <- c(0.12, 0.28, 0.45, 0.64, 0.80, 0.94)
z <- as.matrix(sample_phase[, ..program_cols]); z <- sweep(z, 1, apply(z, 1, max), "-") / 0.8
pp <- exp(z); pp <- pp / rowSums(pp); sample_phase[, phase := as.numeric(pp %*% centers)]
wide <- dcast(sample_phase[n_cells >= 20], donor + coarse_cell_type ~ condition,
              value.var = "phase")
wide[, phase_shift_post_minus_pre := post_3h - pre]
human_sync <- wide[is.finite(phase_shift_post_minus_pre), .(
  response_RSI = regenerative_synchrony_index(phase_shift_post_minus_pre,
                                               sigma_reference = 0.10),
  phase_shift_dispersion = var(phase_shift_post_minus_pre),
  median_phase_shift = median(phase_shift_post_minus_pre),
  effective_cell_types = .N), by = donor]

fwrite(rbindlist(qc_rows), file.path(outdir, "qc_summary.tsv"), sep = "\t")
fwrite(rbindlist(threshold_rows), file.path(outdir, "qc_thresholds.tsv"), sep = "\t")
fwrite(md, file.path(outdir, "cell_metadata.tsv.gz"), sep = "\t")
fwrite(sample_phase, file.path(outdir, "sample_celltype_phase.tsv"), sep = "\t")
fwrite(wide, file.path(outdir, "paired_phase_shifts.tsv"), sep = "\t")
fwrite(human_sync, file.path(outdir, "paired_response_RSI.tsv"), sep = "\t")
metadata(sce)$duplicate_handling <- "GSM6611299/1300 excluded because their barcodes are fully contained in deep-sequenced GSM6611301/1302"
if (!dir.exists(object_out)) {
  HDF5Array::saveHDF5SummarizedExperiment(sce, dir = object_out, replace = FALSE)
} else {
  message("Preserving existing HDF5 checkpoint: ", object_out)
}

p <- ggplot(wide, aes(coarse_cell_type, phase_shift_post_minus_pre, colour = donor,
                      group = donor)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line() + geom_point() + labs(x = NULL, y = "3-h post-exercise phase shift") +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(outdir, "paired_human_phase_shifts.pdf"), p, width = 6.5, height = 3.8)
cat("GSE214544 complete:", ncol(sce), "unique QC singlets across 3 paired donors.\n")

  },
  "spatial" = function() {
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

input_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "extracted/GSE161318")
cache_dir <- file.path(Sys.getenv("MUSCLE_DATA_ROOT", "data"), "metadata/visium_v1_reference")
outdir <- "results/GSE161318_spatial"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Array row/column is a fixed property of each Visium v1/v2 barcode. Cache the
# official 10x example output and use it only as a barcode-to-array lookup.
ref_tar <- file.path(cache_dir, "V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz")
ref_pos <- file.path(cache_dir, "spatial", "tissue_positions_list.csv")
if (!file.exists(ref_pos)) {
  if (!file.exists(ref_tar)) {
    download.file(
      "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz",
      ref_tar, mode = "wb", quiet = FALSE
    )
  }
  untar(ref_tar, exdir = cache_dir)
}
stopifnot(file.exists(ref_pos))
array_map <- fread(ref_pos, header = FALSE,
  col.names = c("barcode", "reference_in_tissue", "array_row", "array_col",
                "reference_pixel_row", "reference_pixel_col"))
stopifnot(nrow(array_map) == 4992L, !anyDuplicated(array_map$barcode))
array_map[, barcode := sub("-1$", "", barcode)]

meta_files <- list.files(input_dir, pattern = "^[^._].*_meta\\.csv\\.gz$", full.names = TRUE)
json_files <- list.files(input_dir, pattern = "^[^._].*\\.json\\.gz$", full.names = TRUE)
stopifnot(length(meta_files) == 3L, length(json_files) == 3L)

stage_columns <- list(
  early_inflammatory = c("theta_M1-Macro.", "theta_Monocyte-(Cxcl10+)",
                         "theta_Monocyte-(Inflammatory)", "theta_Neutrophils",
                         "theta_Dendritic"),
  myogenic_activation = c("theta_Activated-MuSCs", "theta_Committed-Myoblasts"),
  fusion = "theta_Fusing-Myocytes",
  matrix_remodeling = c("theta_FAPs-(Adipogenic)", "theta_FAPs-(Pro-remodeling)",
                        "theta_FAPs-(Stem)", "theta_Tenocytes"),
  mature_repair = c("theta_Myonuclei-(Type-IIb)", "theta_Myonuclei-(Type-IIx)",
                    "theta_Quiescent-MuSCs", "theta_Endothelial-(Capillary)")
)
stage_centers <- c(early_inflammatory = 0.10, myogenic_activation = 0.32,
                   fusion = 0.55, matrix_remodeling = 0.76,
                   mature_repair = 0.95)

read_section <- function(meta_file) {
  x <- fread(meta_file)
  setnames(x, 1L, "barcode")
  x[, barcode_clean := sub("-1$", "", barcode)]
  gsm <- sub("_.*", "", basename(meta_file))
  jf <- json_files[startsWith(basename(json_files), gsm)]
  stopifnot(length(jf) == 1L)
  con <- gzfile(jf, open = "rb"); on.exit(close(con), add = TRUE)
  layout <- fromJSON(con)$oligo
  layout <- as.data.table(layout)[, .(array_row = as.integer(row),
                                     array_col = as.integer(col),
                                     image_x = as.numeric(imageX),
                                     image_y = as.numeric(imageY))]
  stopifnot(nrow(layout) == 4992L,
            !anyDuplicated(layout, by = c("array_row", "array_col")))
  x <- merge(x, array_map[, .(barcode_clean = barcode, array_row, array_col)],
             by = "barcode_clean", all.x = TRUE, sort = FALSE)
  x <- merge(x, layout, by = c("array_row", "array_col"), all.x = TRUE, sort = FALSE)
  if (anyNA(x$image_x) || anyNA(x$image_y)) {
    stop("Incomplete barcode-to-image coordinate mapping for ", gsm)
  }
  x[, section := gsm]
  x
}

spots <- rbindlist(lapply(meta_files, read_section), fill = TRUE)
stopifnot(!anyDuplicated(spots[, paste(section, barcode_clean, sep = "::")]))

# Retain the authors' deconvolution scale: each stage is a sum of mutually
# interpretable theta components and then normalized across the five stages.
for (nm in names(stage_columns)) {
  hit <- intersect(stage_columns[[nm]], names(spots))
  stopifnot(length(hit) == length(stage_columns[[nm]]))
  spots[, (paste0("stage_", nm)) := rowSums(.SD), .SDcols = hit]
}
stage_cols <- paste0("stage_", names(stage_columns))
P <- as.matrix(spots[, ..stage_cols])
P[!is.finite(P) | P < 0] <- 0
P <- P / pmax(rowSums(P), .Machine$double.eps)
spots[, dominant_stage := names(stage_columns)[max.col(P, ties.method = "first")]]
spots[, spatial_phase := as.numeric(P %*% stage_centers)]
spots[, phase_entropy := -rowSums(P * log(pmax(P, 1e-12))) / log(ncol(P))]
spots[, dominant_stage := factor(dominant_stage, levels = names(stage_columns))]

# Hex-grid neighbors use (row +/- 1, col +/- 1) and (row, col +/- 2).
make_edges <- function(d) {
  key <- paste(d$array_row, d$array_col, sep = ":")
  lookup <- setNames(seq_len(nrow(d)), key)
  offsets <- rbind(c(0, 2), c(1, 1), c(1, -1))
  edges <- lapply(seq_len(nrow(d)), function(i) {
    targets <- paste(d$array_row[i] + offsets[, 1], d$array_col[i] + offsets[, 2], sep = ":")
    j <- unname(lookup[targets]); j <- j[!is.na(j)]
    if (!length(j)) return(NULL)
    data.table(i = i, j = j)
  })
  rbindlist(edges)
}

fragmentation <- spots[, {
  d <- .SD
  e <- make_edges(d)
  xi <- d$spatial_phase[e$i]; xj <- d$spatial_phase[e$j]
  xbar <- mean(d$spatial_phase)
  moran <- if (nrow(e) && sum((d$spatial_phase - xbar)^2) > 0) {
    nrow(d) / (2 * nrow(e)) *
      sum(2 * (xi - xbar) * (xj - xbar)) / sum((d$spatial_phase - xbar)^2)
  } else NA_real_
  stage_i <- as.integer(d$dominant_stage[e$i])
  stage_j <- as.integer(d$dominant_stage[e$j])
  .(n_spots = nrow(d), n_edges = nrow(e), mean_phase = mean(d$spatial_phase),
    phase_sd = sd(d$spatial_phase), mean_entropy = mean(d$phase_entropy),
    Moran_I = moran,
    adjacent_stage_boundary_fraction = mean(stage_i != stage_j),
    long_range_phase_boundary_fraction = mean(abs(stage_i - stage_j) > 1))
}, by = .(section, injury)]

fwrite(spots, file.path(outdir, "spatial_phase_spots.tsv.gz"), sep = "\t")
fwrite(fragmentation, file.path(outdir, "spatial_fragmentation_summary.tsv"), sep = "\t")
fwrite(data.table(
  source = c("GSE161318 author spot deconvolution", "10x official Visium v1 reference"),
  use = c("cell-type theta and section-specific image transform",
          "fixed barcode-to-array-row/column lookup only"),
  url = c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161318",
          "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz")
), file.path(outdir, "coordinate_provenance.tsv"), sep = "\t")

stage_palette <- c(early_inflammatory = "#B2182B", myogenic_activation = "#EF8A62",
                   fusion = "#FDB863", matrix_remodeling = "#67A9CF",
                   mature_repair = "#2166AC")
p <- ggplot(spots, aes(image_x, -image_y, colour = dominant_stage)) +
  geom_point(size = 0.42, alpha = 0.9) + facet_wrap(~ injury) +
  scale_colour_manual(values = stage_palette, drop = FALSE) + coord_equal() +
  labs(x = NULL, y = NULL, colour = "Spatial phase") + theme_void(base_size = 9) +
  theme(legend.position = "bottom")
ggsave(file.path(outdir, "spatial_regeneration_phase_map.pdf"), p,
       width = 9.2, height = 3.8)

cat("GSE161318 spatial analysis complete:", nrow(spots), "spots across",
    nrow(fragmentation), "sections.\n")

  },
  "dynamics" = function() {
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(ggplot2)
})

outdir <- "results/complete_dynamics"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
score_cols <- paste0("state_score_", c("Activation", "Proliferation", "Commitment",
  "Fusion", "Remodeling", "Resolution"))
phase_centers <- c(0.12, 0.28, 0.45, 0.64, 0.80, 0.94)

wm <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]; o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}
scale_safe <- function(x) {
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-8) return(rep(0, length(x)))
  (x - m) / s
}

core <- fread("results/discovery_clock/sample_celltype_phase.tsv")
core <- core[, c("dataset_id", "sample_id", "biological_group", "time_days",
  "age_class", "injury_model_std", "condition", "cell_type", "n_cells",
  score_cols), with = FALSE]

read_external <- function(path, dataset, age, injury) {
  x <- fread(path)
  setnames(x, "coarse_cell_type", "cell_type")
  x[, `:=`(dataset_id = dataset, biological_group = sample_id,
           age_class = age, injury_model_std = injury, condition = "")]
  x[, c("dataset_id", "sample_id", "biological_group", "time_days",
        "age_class", "injury_model_std", "condition", "cell_type", "n_cells",
        score_cols), with = FALSE]
}
extra <- rbindlist(list(
  read_external("results/GSE159500_extension/sample_celltype_states.tsv",
                "GSE159500", "young", "NTX"),
  read_external("results/GSE162172_old/sample_celltype_states.tsv",
                "GSE162172", "old", "NTX"),
  read_external("results/GSE205690_bacl2/sample_celltype_states.tsv",
                "GSE205690", "young", "BaCl2")
), fill = TRUE)
x <- rbindlist(list(core, extra), fill = TRUE)
x <- x[n_cells >= 20 & is.finite(time_days)]
x[, sample_uid := paste(dataset_id, sample_id, sep = "::")]
x[, reliability := pmin(1, 1 - exp(-n_cells / 50))]

# Harmonize assay scales without mixing cells: z scores are learned within each
# study from sample-by-lineage pseudobulk summaries.
for (v in score_cols) x[, (paste0("z_", v)) := scale_safe(get(v)), by = dataset_id]
zcols <- paste0("z_", score_cols)
Z <- as.matrix(x[, ..zcols]); Z <- sweep(Z, 1, apply(Z, 1, max), "-") / 0.8
P <- exp(Z); P <- P / rowSums(P)
x[, phase_harmonized := as.numeric(P %*% phase_centers)]
x[, target_tau := pmin(1, log1p(time_days) / log(8))]

sample_tab <- x[, c(
  list(time_days = unique(time_days)[1], target_tau = unique(target_tau)[1],
       age_class = unique(age_class)[1], injury_model_std = unique(injury_model_std)[1],
       biological_group = unique(biological_group)[1],
       consensus_phase = wm(phase_harmonized, reliability),
       phase_locking = Mod(weighted.mean(exp(1i * pi * phase_harmonized), reliability)),
       n_lineages = .N),
  lapply(.SD, function(v) weighted.mean(v, reliability, na.rm = TRUE))),
  by = .(dataset_id, sample_id, sample_uid), .SDcols = zcols]
phase_wide <- dcast(x, dataset_id + sample_id + sample_uid ~ cell_type,
                    value.var = "phase_harmonized", fun.aggregate = median)
sample_tab <- merge(sample_tab, phase_wide, by = c("dataset_id", "sample_id", "sample_uid"),
                    all.x = TRUE)
lineages <- intersect(c("MuSC_myogenic", "macrophage", "FAP", "endothelial",
                        "myonucleus", "neutrophil", "T_NK"), names(sample_tab))
lag_features <- character()
for (i in seq_len(length(lineages) - 1L)) for (j in (i + 1L):length(lineages)) {
  nm <- paste0("lag__", lineages[i], "__", lineages[j])
  sample_tab[, (nm) := sin(pi * (get(lineages[j]) - get(lineages[i])))]
  lag_features <- c(lag_features, nm)
}

ridge_fit <- function(train, features, target = "target_tau", lambda = 1) {
  X <- as.matrix(train[, ..features]); med <- apply(X, 2, median, na.rm = TRUE)
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(X))) X[!is.finite(X[, j]), j] <- med[j]
  mu <- colMeans(X); sc <- apply(X, 2, sd); sc[!is.finite(sc) | sc < 1e-8] <- 1
  X <- sweep(sweep(X, 2, mu, "-"), 2, sc, "/")
  D <- diag(c(0, rep(lambda, ncol(X))))
  if (!identical(dim(crossprod(cbind(1, X))), dim(D))) {
    stop("ridge dimension mismatch for features: ", paste(features, collapse = ","),
         "; X=", paste(dim(X), collapse = "x"),
         "; penalty=", paste(dim(D), collapse = "x"))
  }
  beta <- solve(crossprod(cbind(1, X)) + D,
                crossprod(cbind(1, X), train[[target]]))
  list(features = features, med = med, mu = mu, sc = sc, beta = beta)
}
ridge_predict <- function(fit, test) {
  features <- fit$features
  X <- as.matrix(test[, ..features])
  for (j in seq_len(ncol(X))) X[!is.finite(X[, j]), j] <- fit$med[j]
  X <- sweep(sweep(X, 2, fit$mu, "-"), 2, fit$sc, "/")
  pmin(1.2, pmax(-0.1, as.numeric(cbind(1, X) %*% fit$beta)))
}
choose_lambda <- function(train, features, group = "sample_uid") {
  grid <- c(0.1, 1, 10, 100); ids <- unique(train[[group]])
  if (length(ids) < 6L) return(10)
  errs <- sapply(grid, function(lam) {
    pr <- numeric(); yy <- numeric()
    for (id in ids) {
      tr <- train[get(group) != id]; te <- train[get(group) == id]
      if (nrow(tr) <= length(features) + 2) next
      pr <- c(pr, ridge_predict(ridge_fit(tr, features, lambda = lam), te));
      yy <- c(yy, te$target_tau)
    }
    mean(abs(pr - yy), na.rm = TRUE)
  })
  if (!any(is.finite(errs))) return(10)
  errs[!is.finite(errs)] <- Inf
  grid[which.min(errs)]
}

program_features <- c("consensus_phase", zcols)
independent_predict <- function(train, test) {
  preds <- matrix(NA_real_, nrow(test), length(lineages), dimnames = list(NULL, lineages))
  for (ct in lineages) {
    tr <- train[is.finite(get(ct))]
    if (nrow(tr) < 8L) next
    fit <- ridge_fit(tr, ct, lambda = choose_lambda(tr, ct))
    preds[, ct] <- ridge_predict(fit, test)
  }
  apply(preds, 1, median, na.rm = TRUE)
}
predict_models <- function(train, test) {
  models <- list(
    absolute_time_model = rep(mean(train$target_tau), nrow(test)),
    cell_state_model = ridge_predict(ridge_fit(train, zcols,
      lambda = choose_lambda(train, zcols)), test),
    pseudotime_model = ridge_predict(ridge_fit(train, "consensus_phase",
      lambda = choose_lambda(train, "consensus_phase")), test),
    independent_celltype_clocks = independent_predict(train, test),
    coupled_regeneration_clock = ridge_predict(ridge_fit(train,
      c(lineages, lag_features), lambda = choose_lambda(train, c(lineages, lag_features))), test)
  )
  rbindlist(lapply(names(models), function(m) data.table(
    model = m, sample_uid = test$sample_uid, target_tau = test$target_tau,
    prediction = models[[m]], absolute_error = abs(models[[m]] - test$target_tau))))
}

reference <- sample_tab[age_class == "young" & injury_model_std %in% c("NTX", "CTX", "BaCl2")]
run_folds <- function(fold_name, fold_vector, data = reference) {
  out <- list(); vals <- unique(fold_vector)
  for (v in vals) {
    te <- which(fold_vector == v); tr <- which(fold_vector != v)
    if (length(te) < 1L || length(tr) < 8L) next
    z <- predict_models(data[tr], data[te]); z[, `:=`(validation = fold_name, held_out = as.character(v))]
    out[[as.character(v)]] <- z
  }
  rbindlist(out, fill = TRUE)
}
cv <- rbindlist(list(
  run_folds("leave_study_out", reference$dataset_id),
  run_folds("leave_time_point_out", as.character(reference$time_days)),
  run_folds("leave_injury_model_out", reference$injury_model_std)
), fill = TRUE)

# Leave-age-out is intentionally directional: fit successful young repair and
# evaluate old/geriatric samples without refitting age-specific coefficients.
age_test <- sample_tab[age_class %in% c("old", "geriatric") & injury_model_std == "NTX"]
if (nrow(age_test)) {
  age_cv <- predict_models(reference[injury_model_std == "NTX"], age_test)
  age_cv[, `:=`(validation = "leave_age_out", held_out = age_test$age_class[match(sample_uid, age_test$sample_uid)])]
  cv <- rbindlist(list(cv, age_cv), fill = TRUE)
}
cv_metrics <- cv[, .(n_samples = sum(is.finite(absolute_error)),
  MAE_tau = mean(absolute_error, na.rm = TRUE),
  RMSE_tau = sqrt(mean((prediction - target_tau)^2, na.rm = TRUE))),
  by = .(validation, held_out, model)]

# Leave-cell-type-out ablation of the coupled clock.
celltype_ablation <- rbindlist(lapply(lineages, function(drop) {
  features <- c(setdiff(lineages, drop), lag_features[!grepl(paste0("(^|__)", drop, "(__|$)"), lag_features)])
  pr <- numeric(nrow(reference))
  for (study in unique(reference$dataset_id)) {
    tr <- reference[dataset_id != study]; te <- reference[dataset_id == study]
    fit <- ridge_fit(tr, features, lambda = choose_lambda(tr, features))
    pr[reference$dataset_id == study] <- ridge_predict(fit, te)
  }
  data.table(removed_cell_type = drop, MAE_tau = mean(abs(pr - reference$target_tau)),
             delta_MAE_vs_full = NA_real_)
}))
full_mae <- cv_metrics[validation == "leave_study_out" & model == "coupled_regeneration_clock",
                       weighted.mean(MAE_tau, n_samples)]
celltype_ablation[, delta_MAE_vs_full := MAE_tau - full_mae]

# Dynamic time warping on multivariate lineage-phase trajectories.
dtw_distance <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B); n <- nrow(A); m <- nrow(B)
  D <- matrix(Inf, n + 1L, m + 1L); D[1, 1] <- 0
  for (i in seq_len(n)) for (j in seq_len(m)) {
    cost <- sqrt(mean((A[i, ] - B[j, ])^2, na.rm = TRUE))
    if (!is.finite(cost)) cost <- 1
    D[i + 1L, j + 1L] <- cost + min(D[i, j + 1L], D[i + 1L, j], D[i, j])
  }
  D[n + 1L, m + 1L] / (n + m)
}
study_curves <- sample_tab[age_class == "young", lapply(.SD, median, na.rm = TRUE),
                           by = .(dataset_id, time_days), .SDcols = lineages]
study_pairs <- combn(unique(study_curves$dataset_id), 2, simplify = FALSE)
dtw <- rbindlist(lapply(study_pairs, function(pp) {
  a <- study_curves[dataset_id == pp[1]][order(time_days), ..lineages]
  b <- study_curves[dataset_id == pp[2]][order(time_days), ..lineages]
  data.table(study_a = pp[1], study_b = pp[2], multivariate_DTW_cost = dtw_distance(a, b),
             n_time_a = nrow(a), n_time_b = nrow(b))
}))

# Exact one-dimensional Wasserstein-1 distance between empirical phase distributions.
w1 <- function(a, b, nq = 101L) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  mean(abs(quantile(a, seq(0, 1, length.out = nq), names = FALSE) -
             quantile(b, seq(0, 1, length.out = nq), names = FALSE)))
}
ot <- sample_tab[age_class == "young", {
  tt <- sort(unique(time_days)); if (length(tt) < 2L) return(NULL)
  rbindlist(lapply(seq_len(length(tt) - 1L), function(k) {
    a <- unlist(.SD[time_days == tt[k], lineages, with = FALSE])
    b <- unlist(.SD[time_days == tt[k + 1L], lineages, with = FALSE])
    data.table(time_from = tt[k], time_to = tt[k + 1L], delta_time = tt[k + 1L] - tt[k],
               Wasserstein_1 = w1(a, b), transport_speed = w1(a, b) / (tt[k + 1L] - tt[k]))
  }))
}, by = dataset_id]

# GP/GAM reference trajectories and timing abnormality calls.
x[, expected_phase := NA_real_]; x[, reference_sd := NA_real_]
young_ref <- x[age_class == "young" & injury_model_std %in% c("NTX", "CTX", "BaCl2")]
gp_landmarks <- list()
for (ct in unique(young_ref$cell_type)) {
  tr <- young_ref[cell_type == ct]
  if (nrow(tr) < 10L || uniqueN(tr$target_tau) < 4L) next
  k <- min(5L, uniqueN(tr$target_tau) - 1L)
  fit <- gam(phase_harmonized ~ s(target_tau, k = k),
             data = tr, weights = reliability, method = "REML")
  idx <- which(x$cell_type == ct)
  x$expected_phase[idx] <- predict(fit, newdata = x[idx])
  rsd <- sd(residuals(fit), na.rm = TRUE); x$reference_sd[idx] <- rsd
  grid <- data.table(target_tau = seq(0, 1, length.out = 201))
  pg <- predict(fit, newdata = grid)
  gp_landmarks[[ct]] <- data.table(cell_type = ct,
    peak_tau = grid$target_tau[which.max(pg)], peak_phase = max(pg),
    exit_tau = grid$target_tau[max(which(pg >= min(pg) + 0.25 * diff(range(pg))))],
    reference_residual_sd = rsd)
}
x[, phase_deviation := phase_harmonized - expected_phase]
x[, sigma_phase_reference := max(0.05,
  sd(phase_deviation[age_class == "young" & injury_model_std %in% c("NTX", "CTX", "BaCl2")],
     na.rm = TRUE))]
x[, timing_class := fifelse(!is.finite(phase_deviation) | !is.finite(reference_sd), NA_character_,
  fifelse(phase_deviation > 1.96 * reference_sd, "phase_advance",
    fifelse(phase_deviation < -1.96 * reference_sd, "phase_delay", "within_reference")))]
abnormal <- x[age_class %in% c("old", "geriatric"), .(
  n_samples = uniqueN(sample_uid), median_deviation = median(phase_deviation, na.rm = TRUE),
  phase_advance_fraction = mean(timing_class == "phase_advance", na.rm = TRUE),
  phase_delay_fraction = mean(timing_class == "phase_delay", na.rm = TRUE)),
  by = .(age_class, time_days, cell_type)]

# Young-vs-aged synchrony and critical-window analysis within NTX.
rsi_age <- x[injury_model_std == "NTX" & is.finite(phase_deviation), .(
  RSI = regenerative_synchrony_index(phase_deviation, reliability,
    sigma_reference = unique(sigma_phase_reference)[1]),
  dispersion = weighted_var(phase_deviation, reliability), effective_cell_types = .N),
  by = .(dataset_id, sample_id, sample_uid, age_class, time_days)]
window <- rsi_age[, .(mean_RSI = mean(RSI, na.rm = TRUE),
  mean_dispersion = mean(dispersion, na.rm = TRUE), n_samples = .N),
  by = .(age_class, time_days)]
window_wide <- dcast(window, time_days ~ age_class, value.var = "mean_RSI")
if (all(c("young", "geriatric") %in% names(window_wide)))
  window_wide[, geriatric_minus_young := geriatric - young]

# Linear Gaussian state update pooled across young studies; evaluation is
# leave-study-out future condition-level prediction.
make_transitions <- function(d) {
  curves <- d[, lapply(.SD, median, na.rm = TRUE), by = .(dataset_id, time_days), .SDcols = lineages]
  rbindlist(lapply(unique(curves$dataset_id), function(ds) {
    z <- curves[dataset_id == ds][order(time_days)]; if (nrow(z) < 2L) return(NULL)
    a <- z[-.N]; b <- z[-1L]
    out <- data.table(dataset_id = ds, time_from = a$time_days, time_to = b$time_days,
                      delta_time = b$time_days - a$time_days)
    from <- copy(a[, ..lineages]); setnames(from, paste0("from_", lineages))
    to <- copy(b[, ..lineages]); setnames(to, paste0("to_", lineages))
    cbind(out, from, to)
  }), fill = TRUE)
}
transitions <- make_transitions(sample_tab[age_class == "young"])
state_cv <- rbindlist(lapply(unique(transitions$dataset_id), function(held) {
  tr <- transitions[dataset_id != held]; te <- transitions[dataset_id == held]
  out <- list()
  for (target in lineages) {
    feats <- c(paste0("from_", lineages), "delta_time")
    tr2 <- tr[is.finite(get(paste0("to_", target)))]
    if (nrow(tr2) < 6L) next
    fit <- ridge_fit(tr2, feats, target = paste0("to_", target), lambda = 10)
    pred <- ridge_predict(fit, te)
    out[[target]] <- data.table(held_out_study = held, cell_type = target,
      time_from = te$time_from, time_to = te$time_to,
      observed = te[[paste0("to_", target)]], predicted = pred,
      squared_error = (pred - te[[paste0("to_", target)]])^2)
  }
  rbindlist(out, fill = TRUE)
}), fill = TRUE)
state_metrics <- state_cv[, .(n_transitions = sum(is.finite(squared_error)),
  future_state_RMSE = sqrt(mean(squared_error, na.rm = TRUE))), by = .(held_out_study, cell_type)]

# Kuramoto-style phenomenological coupling, without wrapping the transient phase.
estimate_K <- function(d, age_label) {
  curve <- d[, lapply(.SD, median, na.rm = TRUE), by = time_days, .SDcols = lineages][order(time_days)]
  for (v in lineages) {
    vv <- curve[[v]]; ok <- is.finite(vv)
    if (sum(ok) >= 2L) vv <- approx(curve$time_days[ok], vv[ok], curve$time_days,
                                    rule = 2, ties = "ordered")$y
    if (any(!is.finite(vv))) vv[!is.finite(vv)] <- median(vv, na.rm = TRUE)
    curve[[v]] <- vv
  }
  curve <- curve[complete.cases(curve[, ..lineages])]
  if (nrow(curve) < 4L) return(NULL)
  theta <- as.matrix(curve[, ..lineages]) * pi; dt <- diff(curve$time_days)
  deriv <- sweep(diff(theta), 1, dt, "/"); mid <- (theta[-1, ] + theta[-nrow(theta), ]) / 2
  rbindlist(lapply(seq_along(lineages), function(i) {
    others <- setdiff(seq_along(lineages), i)
    Xk <- cbind(1, sapply(others, function(j) sin(mid[, j] - mid[, i])))
    beta <- solve(crossprod(Xk) + diag(c(0, rep(1, length(others)))), crossprod(Xk, deriv[, i]))
    data.table(age_class = age_label, target = lineages[i], source = lineages[others], K = beta[-1])
  }))
}
K <- rbindlist(list(
  estimate_K(sample_tab[age_class == "young" & injury_model_std == "NTX"], "young"),
  estimate_K(sample_tab[age_class %in% c("old", "geriatric") & injury_model_std == "NTX"], "aged")
), fill = TRUE)

# Human leave-species transfer: the clock is trained only in mouse and tested
# on paired pre/post human samples. Outcome is ordering, not toxin-repair MAE.
human <- fread("results/GSE214544_human/sample_celltype_phase.tsv")
human_order <- human[, .(consensus_phase = median(phase, na.rm = TRUE), n_lineages = .N),
                     by = .(donor, condition)]
human_wide <- dcast(human_order, donor ~ condition, value.var = "consensus_phase")
human_wide[, post_minus_pre := post_3h - pre]
species_metric <- data.table(validation = "leave_species_out",
  n_donors = nrow(human_wide), directional_concordance = mean(human_wide$post_minus_pre > 0),
  median_phase_shift = median(human_wide$post_minus_pre))

fwrite(x, file.path(outdir, "harmonized_sample_celltype_atlas.tsv.gz"), sep = "\t")
fwrite(sample_tab, file.path(outdir, "harmonized_sample_clock.tsv"), sep = "\t")
fwrite(cv, file.path(outdir, "all_validation_predictions.tsv.gz"), sep = "\t")
fwrite(cv_metrics, file.path(outdir, "all_validation_metrics.tsv"), sep = "\t")
fwrite(celltype_ablation, file.path(outdir, "leave_cell_type_out.tsv"), sep = "\t")
fwrite(dtw, file.path(outdir, "DTW_alignment_costs.tsv"), sep = "\t")
fwrite(ot, file.path(outdir, "optimal_transport_transitions.tsv"), sep = "\t")
fwrite(rbindlist(gp_landmarks), file.path(outdir, "GP_phase_landmarks.tsv"), sep = "\t")
fwrite(abnormal, file.path(outdir, "age_abnormal_timing_phenotypes.tsv"), sep = "\t")
fwrite(rsi_age, file.path(outdir, "age_RSI_all_samples.tsv"), sep = "\t")
fwrite(window_wide, file.path(outdir, "critical_window_age_effect.tsv"), sep = "\t")
fwrite(state_cv, file.path(outdir, "state_space_future_predictions.tsv"), sep = "\t")
fwrite(state_metrics, file.path(outdir, "state_space_future_metrics.tsv"), sep = "\t")
fwrite(K, file.path(outdir, "Kuramoto_coupling_by_age.tsv"), sep = "\t")
fwrite(human_wide, file.path(outdir, "human_leave_species_out.tsv"), sep = "\t")
fwrite(species_metric, file.path(outdir, "human_leave_species_metric.tsv"), sep = "\t")

p <- ggplot(cv_metrics[validation != "leave_age_out"],
            aes(model, MAE_tau, colour = validation)) + geom_point(position = position_jitter(width = .12)) +
  coord_flip() + theme_classic(base_size = 9) + labs(x = NULL, y = "Biological-time MAE")
ggsave(file.path(outdir, "complete_validation_benchmark.pdf"), p, width = 7.2, height = 4.3)
p2 <- ggplot(rsi_age[is.finite(RSI)], aes(time_days, RSI, colour = age_class)) + geom_point(alpha = .65) +
  stat_summary(fun = mean, geom = "line", aes(group = age_class), linewidth = .8) +
  theme_classic(base_size = 9) + labs(x = "Days post injury", y = "RSI")
ggsave(file.path(outdir, "age_RSI_complete_timecourse.pdf"), p2, width = 6.5, height = 3.8)

writeLines(c(
  "DTW: exact dynamic-programming multivariate alignment on sample-level lineage phases",
  "OT: exact empirical one-dimensional Wasserstein-1 distance between adjacent phase distributions",
  "GP: mgcv REML smooths with study fixed effects and sample-lineage reliability weights",
  "state_space: ridge-regularized linear Gaussian mean transition; condition-level future prediction",
  "Kuramoto: unwrapped transient phenomenological coupling; coefficients are predictive, not causal",
  "all validation folds keep complete biological samples together"
), file.path(outdir, "method_manifest.txt"))
cat("Complete dynamical validation finished:", nrow(sample_tab), "samples,",
    nrow(cv), "held-out predictions.\n")

  },
  "programmes" = function() {
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(ggplot2)
})

indir <- "results/complete_dynamics"
outdir <- "results/program_timing"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(indir, "harmonized_sample_celltype_atlas.tsv.gz"))
programs <- c("Activation", "Proliferation", "Commitment", "Fusion",
              "Remodeling", "Resolution")
zcols <- paste0("z_state_score_", programs)
x <- x[injury_model_std == "NTX" & age_class %in% c("young", "old", "geriatric")]
x[, age_binary := as.integer(age_class != "young")]
x[, age_group2 := ifelse(age_binary == 1L, "aged", "young")]

decomp <- list(); curve_rows <- list()
for (ct in unique(x$cell_type)) for (ii in seq_along(programs)) {
  ycol <- zcols[ii]; d <- x[cell_type == ct & is.finite(get(ycol))]
  if (nrow(d) < 18L || uniqueN(d$target_tau) < 4L || uniqueN(d$age_binary) < 2L) next
  d[, y := get(ycol)]
  degree <- min(3L, d[, uniqueN(target_tau), by = age_binary]$V1 - 1L)
  if (degree < 1L) next
  null <- try(lm(y ~ age_binary + poly(target_tau, degree = degree, raw = TRUE),
                 data = d, weights = reliability), silent = TRUE)
  full <- try(lm(y ~ age_binary * poly(target_tau, degree = degree, raw = TRUE),
                 data = d, weights = reliability), silent = TRUE)
  if (inherits(null, "try-error") || inherits(full, "try-error")) next
  aa <- try(anova(null, full), silent = TRUE)
  timing_p <- if (!inherits(aa, "try-error") && nrow(aa) >= 2L) aa$`Pr(>F)`[2] else NA_real_
  grid <- rbindlist(list(
    data.table(target_tau = seq(0, 1, length.out = 201), age_binary = 0L),
    data.table(target_tau = seq(0, 1, length.out = 201), age_binary = 1L)
  ))
  grid[, predicted := as.numeric(predict(full, newdata = grid))]
  py <- grid[age_binary == 0L]; pa <- grid[age_binary == 1L]
  peak_y <- py$target_tau[which.max(py$predicted)]
  peak_a <- pa$target_tau[which.max(pa$predicted)]
  peak_amp <- max(pa$predicted) - max(py$predicted)
  aligned_amp <- mean(pa$predicted - py$predicted)
  peak_shift <- peak_a - peak_y
  cls <- if (abs(peak_shift) >= 0.10 && abs(peak_amp) < 0.25) "timing_dominant" else
    if (abs(peak_shift) < 0.10 && abs(peak_amp) >= 0.25) "strength_dominant" else
      if (abs(peak_shift) >= 0.10 && abs(peak_amp) >= 0.25) "mixed" else "conserved"
  key <- paste(ct, programs[ii], sep = "::")
  decomp[[key]] <- data.table(cell_type = ct, program = programs[ii], n_rows = nrow(d),
    aged_aligned_amplitude = aligned_amp, aged_peak_amplitude_difference = peak_amp,
    aged_peak_shift_tau = peak_shift, timing_interaction_p = timing_p,
    dominant_abnormality = cls)
  grid[, `:=`(cell_type = ct, program = programs[ii],
              age_group = ifelse(age_binary == 1L, "aged", "young"))]
  curve_rows[[key]] <- grid
}
decomp <- rbindlist(decomp, fill = TRUE)
curves <- rbindlist(curve_rows, fill = TRUE)

# Frozen young-reference envelopes for prespecified late persistence phenotypes.
tests <- data.table(
  cell_type = c("macrophage", "macrophage", "FAP", "MuSC_myogenic", "MuSC_myogenic"),
  program = c("Activation", "Resolution", "Remodeling", "Activation", "Fusion"),
  phenotype = c("inflammatory_persistence", "failed_resolution", "ECM_persistence",
                "ectopic_activation", "delayed_differentiation")
)
persistence <- list()
for (i in seq_len(nrow(tests))) {
  tt <- tests[i]; ycol <- paste0("z_state_score_", tt$program)
  ref <- x[cell_type == tt$cell_type & age_binary == 0L & is.finite(get(ycol))]
  aged <- x[cell_type == tt$cell_type & age_binary == 1L & time_days >= 5 & is.finite(get(ycol))]
  if (nrow(ref) < 8L || !nrow(aged)) next
  ref[, y := get(ycol)]; k <- min(4L, uniqueN(ref$target_tau) - 1L)
  if (k < 2L) next
  fit <- lm(y ~ poly(target_tau, degree = min(2L, k), raw = TRUE),
            data = ref, weights = reliability)
  aged[, expected := predict(fit, newdata = aged)]
  rsd <- sd(residuals(fit), na.rm = TRUE)
  aged[, abnormal := if (tt$phenotype == "delayed_differentiation")
    get(ycol) < expected - 1.96 * rsd else get(ycol) > expected + 1.96 * rsd]
  persistence[[i]] <- aged[, .(cell_type = tt$cell_type, program = tt$program,
    phenotype = tt$phenotype, time_days, age_class, sample_uid,
    observed = get(ycol), expected, reference_sd = rsd, abnormal)]
}
persistence <- rbindlist(persistence, fill = TRUE)
persistence_summary <- persistence[, .(n_samples = .N,
  abnormal_fraction = mean(abnormal, na.rm = TRUE),
  median_excess = median(observed - expected, na.rm = TRUE)),
  by = .(phenotype, cell_type, program, time_days, age_class)]

# Phase-lag curves and temporal compression across the three locked relationships.
s <- fread(file.path(indir, "harmonized_sample_clock.tsv"))
pairs <- data.table(a = c("MuSC_myogenic", "FAP", "endothelial"),
                    b = c("macrophage", "MuSC_myogenic", "myonucleus"),
                    contrast = c("MuSC_minus_macrophage", "FAP_minus_MuSC",
                                 "endothelial_minus_myonucleus"))
lag_rows <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  pp <- pairs[i]
  if (!all(c(pp$a, pp$b) %in% names(s))) return(NULL)
  s[, .(dataset_id, sample_uid, age_class, injury_model_std, time_days,
        contrast = pp$contrast, phase_lag = get(pp$a) - get(pp$b))]
}))
lag_summary <- lag_rows[injury_model_std == "NTX" & age_class %in% c("young", "old", "geriatric"),
  .(mean_lag = mean(phase_lag, na.rm = TRUE), sd_lag = sd(phase_lag, na.rm = TRUE),
    n_samples = sum(is.finite(phase_lag))), by = .(contrast, age_class, time_days)]

# Peak spacing is a transparent temporal-compression descriptor.
compression <- data.table(
  metric = c("aged_activation_to_fusion_peak_shift_difference",
             "aged_commitment_to_fusion_peak_shift_difference"),
  value = c(decomp[cell_type == "MuSC_myogenic" & program == "Fusion", aged_peak_shift_tau] -
              decomp[cell_type == "MuSC_myogenic" & program == "Activation", aged_peak_shift_tau],
            decomp[cell_type == "MuSC_myogenic" & program == "Fusion", aged_peak_shift_tau] -
              decomp[cell_type == "MuSC_myogenic" & program == "Commitment", aged_peak_shift_tau])
)

fwrite(decomp, file.path(outdir, "program_strength_vs_timing.tsv"), sep = "\t")
fwrite(curves, file.path(outdir, "fitted_program_curves.tsv.gz"), sep = "\t")
fwrite(persistence, file.path(outdir, "late_abnormal_program_calls.tsv"), sep = "\t")
fwrite(persistence_summary, file.path(outdir, "late_abnormal_program_summary.tsv"), sep = "\t")
fwrite(lag_rows, file.path(outdir, "prespecified_phase_lags_all_samples.tsv"), sep = "\t")
fwrite(lag_summary, file.path(outdir, "prespecified_phase_lag_summary.tsv"), sep = "\t")
fwrite(compression, file.path(outdir, "temporal_compression_metrics.tsv"), sep = "\t")

p <- ggplot(decomp, aes(aged_peak_shift_tau, aged_peak_amplitude_difference,
                        colour = dominant_abnormality)) +
  geom_hline(yintercept = c(-.25, .25), colour = "grey80", linewidth = .3) +
  geom_vline(xintercept = c(-.1, .1), colour = "grey80", linewidth = .3) +
  geom_point() + facet_wrap(~cell_type) + theme_classic(base_size = 8) +
  labs(x = "Aged peak shift in biological time", y = "Aged peak-amplitude difference")
ggsave(file.path(outdir, "program_timing_amplitude_decomposition.pdf"), p,
       width = 9, height = 6)
cat("Program timing/strength decomposition complete:", nrow(decomp), "lineage-program models.\n")

  },
  "curio" = function() {
suppressPackageStartupMessages({
  library(zellkonverter)
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

indir <- "data/GSE266933"
outdir <- "results/GSE266933_age_spatial_expression"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(indir, pattern = "^[^._].*\\.h5ad$", full.names = TRUE)
stopifnot(length(files) == 2L)

programs <- list(
  AP1 = c("Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1"),
  NFkB = c("Rela", "Relb", "Nfkb1", "Nfkb2", "Nfkbia", "Tnfaip3"),
  interferon = c("Isg15", "Ifit1", "Ifit2", "Ifit3", "Irf7", "Stat1"),
  TGFb = c("Tgfb1", "Tgfbr1", "Smad2", "Smad3", "Ctgf", "Serpine1"),
  STAT = c("Stat1", "Stat3", "Stat5a", "Socs1", "Socs3"),
  MYOD_MYOG = c("Myod1", "Myog", "Myf5", "Myf6"),
  Notch = c("Notch1", "Notch2", "Hes1", "Hey1", "Dll1", "Jag1"),
  Wnt = c("Ctnnb1", "Axin2", "Lef1", "Wnt4", "Wnt7a", "Fzd7"),
  ECM = c("Col1a1", "Col1a2", "Col3a1", "Fn1", "Dcn", "Lum", "Mmp2", "Mmp14"),
  angiogenesis = c("Pecam1", "Kdr", "Cdh5", "Emcn", "Klf2", "Vegfa", "Angpt1"),
  senescence = c("Cdkn1a", "Cdkn2a", "Trp53", "Glb1", "Serpine1", "Gdf15"),
  fusion = c("Myog", "Mymk", "Mymx", "Tnnt1", "Ckm")
)

read_score <- function(f) {
  sce <- readH5AD(f, reader = "R")
  m <- assay(sce, "X")
  age <- if (grepl("Geriatric", basename(f), fixed = TRUE)) "geriatric" else "young"
  totals <- Matrix::colSums(m); detected <- Matrix::colSums(m > 0)
  pos <- totals > 0
  low_total <- max(5, as.numeric(quantile(totals[pos], 0.10, na.rm = TRUE)))
  low_detected <- max(3, as.numeric(quantile(detected[pos], 0.10, na.rm = TRUE)))
  keep <- totals >= low_total & detected >= low_detected
  dt <- data.table(bead_id = colnames(sce), age_class = age,
                   total_counts = totals, detected_genes = detected,
                   passes_expression_qc = keep)
  for (nm in names(programs)) {
    hit <- intersect(programs[[nm]], rownames(sce))
    stopifnot(length(hit) >= 2L)
    dt[[paste0("program_", nm)]] <- log1p(Matrix::colSums(m[hit, , drop = FALSE]) /
                                            pmax(totals, 1) * 1e4 / length(hit))
  }
  qc <- data.table(age_class = age, input_beads = ncol(sce),
                   positive_beads = sum(pos), retained_beads = sum(keep),
                   low_total = low_total, low_detected = low_detected,
                   coordinate_fields = 0L)
  rm(sce, m); gc()
  list(metadata = dt, qc = qc)
}

z <- lapply(files, read_score)
beads <- rbindlist(lapply(z, `[[`, "metadata"))
qc <- rbindlist(lapply(z, `[[`, "qc"))
beads <- beads[passes_expression_qc == TRUE]
pcols <- paste0("program_", names(programs))
for (v in pcols) {
  mu <- mean(beads[[v]], na.rm = TRUE); ss <- sd(beads[[v]], na.rm = TRUE)
  beads[[paste0("z_", v)]] <- (beads[[v]] - mu) / ifelse(ss > 0, ss, 1)
}

stage_defs <- list(
  early_inflammatory = c("z_program_AP1", "z_program_NFkB", "z_program_interferon", "z_program_STAT"),
  myogenic_activation = c("z_program_MYOD_MYOG", "z_program_Notch"),
  fusion = "z_program_fusion",
  matrix_remodeling = c("z_program_TGFb", "z_program_ECM"),
  mature_repair = c("z_program_angiogenesis", "z_program_Wnt")
)
stage_names <- names(stage_defs)
S <- sapply(stage_defs, function(v) rowMeans(as.matrix(beads[, ..v])))
S <- sweep(S, 1, apply(S, 1, max), "-") / 0.8
P <- exp(S); P <- P / rowSums(P)
centers <- c(.10, .32, .55, .76, .95)
set(beads, j = "dominant_stage", value = stage_names[max.col(P, ties.method = "first")])
set(beads, j = "phase", value = as.numeric(P %*% centers))
set(beads, j = "phase_entropy",
    value = -rowSums(P * log(pmax(P, 1e-12))) / log(ncol(P)))

phase_summary <- beads[, .(n_beads = .N, mean_phase = mean(phase),
  phase_sd = sd(phase), mean_entropy = mean(phase_entropy),
  q10_phase = quantile(phase, .1), median_phase = median(phase), q90_phase = quantile(phase, .9)),
  by = age_class]
program_summary <- beads[, lapply(.SD, mean), by = age_class, .SDcols = pcols]
program_long <- melt(program_summary, id.vars = "age_class", variable.name = "program",
                     value.name = "mean_log_normalized_activity")
program_wide <- dcast(program_long, program ~ age_class,
                      value.var = "mean_log_normalized_activity")
program_wide[, geriatric_minus_young := geriatric - young]

stage_composition <- beads[, .N, by = .(age_class, dominant_stage)]
stage_composition[, fraction := N / sum(N), by = age_class]

fwrite(qc, file.path(outdir, "bead_QC_summary.tsv"), sep = "\t")
fwrite(beads, file.path(outdir, "bead_program_phase_scores.tsv.gz"), sep = "\t")
fwrite(phase_summary, file.path(outdir, "age_phase_distribution_summary.tsv"), sep = "\t")
fwrite(program_wide, file.path(outdir, "age_program_effect_sizes.tsv"), sep = "\t")
fwrite(stage_composition, file.path(outdir, "age_stage_composition.tsv"), sep = "\t")
writeLines(c(
  "platform: Curio Seeker spatial beads, one young and one geriatric D5 section",
  "coordinates: absent from both GEO H5AD obs/obsm; no spatial adjacency or patch claims",
  "QC: expression-bearing beads filtered by within-section 10th percentiles, with conservative floors",
  "doublets: scDblFinder not applicable to subcellular spatial bead measurements",
  "inference: descriptive section-level effect sizes only; beads are not biological replicates",
  "phase_fragmentation: distributional heterogeneity only, not spatial fragmentation"
), file.path(outdir, "analysis_manifest.txt"))

p <- ggplot(beads, aes(phase, colour = age_class, fill = age_class)) +
  geom_density(alpha = .16) + theme_classic(base_size = 9) +
  labs(x = "Expression-only regeneration phase", y = "Density")
ggsave(file.path(outdir, "young_geriatric_phase_distribution.pdf"), p,
       width = 5.8, height = 3.7)
cat("GSE266933 expression-only age validation complete:", nrow(beads), "retained beads.\n")

  },
  "human_age" = function() {
suppressPackageStartupMessages({
  library(DropletUtils)
  library(SingleCellExperiment)
  library(scuttle)
  library(scDblFinder)
  library(BiocParallel)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

indir <- "data/GSE167186"
outdir <- "results/GSE167186_human_age"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(indir, pattern = "^[^._].*[.]h5$", full.names = TRUE)
stopifnot(length(files) == 17L)
young_hm <- c(5, 8, 10, 13, 21, 22)

read_one <- function(f) {
  sce <- read10xCounts(f, col.names = TRUE)
  sym <- rowData(sce)$Symbol; sym[is.na(sym) | sym == ""] <- rownames(sce)[is.na(sym) | sym == ""]
  rownames(sce) <- make.unique(sym)
  sid <- sub("_filtered.*", "", basename(f))
  hm <- as.integer(sub(".*_HM([0-9]+).*", "\\1", basename(f)))
  sce$sample_id <- sid; sce$donor <- paste0("HM", hm)
  sce$age_class <- ifelse(hm %in% young_hm, "young", "old")
  colnames(sce) <- paste(sid, colnames(sce), sep = "::")
  sce
}

set.seed(20260825)
objects <- lapply(files, read_one)
qc_rows <- list(); threshold_rows <- list(); all_qc_rows <- list()
for (i in seq_along(objects)) {
  sce <- objects[[i]]; sid <- unique(sce$sample_id)
  message("QC/doublet: ", sid, " (", ncol(sce), " cells)")
  qc <- perCellQCMetrics(sce, subsets = list(Mt = grepl("^MT-", rownames(sce))))
  robust_low <- function(v) max(1, median(v) - 3 * mad(v))
  robust_high <- function(v) median(v) + 5 * mad(v)
  low_detected <- robust_low(qc$detected); low_sum <- robust_low(qc$sum)
  high_detected <- robust_high(qc$detected); high_sum <- robust_high(qc$sum)
  high_mt <- min(40, max(median(qc$subsets_Mt_percent) + 3 * mad(qc$subsets_Mt_percent),
                         quantile(qc$subsets_Mt_percent, 0.99, na.rm = TRUE)))
  pass <- qc$detected >= low_detected & qc$sum >= low_sum &
    qc$detected <= high_detected & qc$sum <= high_sum & qc$subsets_Mt_percent <= high_mt
  if (sum(pass) < 50L) stop("Too few QC-pass cells in ", sid, ": ", sum(pass))
  colData(sce) <- cbind(colData(sce), qc); sce$passes_core_qc <- pass
  dbl <- scDblFinder(sce[, pass], samples = rep(sid, sum(pass)),
                     BPPARAM = SerialParam(), verbose = FALSE)
  cls <- rep(NA_character_, ncol(sce)); score <- rep(NA_real_, ncol(sce))
  cls[pass] <- as.character(dbl$scDblFinder.class); score[pass] <- dbl$scDblFinder.score
  sce$scDblFinder.class <- cls; sce$scDblFinder.score <- score
  sce$passes_qc_and_singlet <- pass & cls == "singlet"
  all_qc_rows[[sid]] <- data.table(cell_id = colnames(sce), sample_id = sid,
    donor = unique(sce$donor), age_class = unique(sce$age_class),
    total_counts = qc$sum, detected_genes = qc$detected,
    mitochondrial_percent = qc$subsets_Mt_percent, passes_core_qc = pass,
    scDblFinder_class = cls, scDblFinder_score = score,
    passes_qc_and_singlet = sce$passes_qc_and_singlet)
  qc_rows[[sid]] <- data.table(sample_id = sid, donor = unique(sce$donor),
    age_class = unique(sce$age_class), input_cells = ncol(sce), core_qc_pass = sum(pass),
    doublets = sum(cls == "doublet", na.rm = TRUE),
    final_singlets = sum(sce$passes_qc_and_singlet, na.rm = TRUE))
  threshold_rows[[sid]] <- data.table(sample_id = sid, low_detected, low_sum,
    high_detected, high_sum, high_mt)
  objects[[i]] <- sce[, sce$passes_qc_and_singlet %in% TRUE]
}

sce <- do.call(cbind, objects); rm(objects); gc()
sce <- logNormCounts(sce)
markers <- list(
  MuSC_myogenic = c("PAX7", "CALCR", "SDC4", "VCAM1", "MYOD1", "MYOG", "DES"),
  myonucleus = c("TTN", "MYH1", "MYH2", "MYH4", "CKM", "ACTA1"),
  FAP = c("PDGFRA", "COL1A1", "COL3A1", "DCN", "LUM"),
  macrophage = c("PTPRC", "LYZ", "CD68", "C1QA", "C1QB"),
  neutrophil = c("S100A8", "S100A9", "FCGR3B", "CSF3R"),
  T_NK = c("CD3D", "CD3E", "TRBC1", "NKG7", "KLRD1"),
  endothelial = c("PECAM1", "KDR", "EMCN", "CDH5", "KLF2"),
  pericyte = c("RGS5", "CSPG4", "PDGFRB", "NOTCH3", "DES")
)
programs <- list(
  Activation = c("FOS", "JUN", "JUNB", "EGR1", "MYOD1"),
  Proliferation = c("MKI67", "TOP2A", "CCNB1", "PCNA", "UBE2C"),
  Commitment = c("MYOD1", "MYF5", "DES", "CDKN1C"),
  Fusion = c("MYOG", "MYMK", "MYMX", "TNNT1", "CKM"),
  Remodeling = c("COL1A1", "COL3A1", "FN1", "MMP14", "TGFB1"),
  Resolution = c("MRC1", "C1QA", "C1QB", "IL10", "APOE")
)
score_module <- function(gs) {
  hit <- intersect(gs, rownames(sce))
  if (!length(hit)) return(rep(NA_real_, ncol(sce)))
  Matrix::colMeans(logcounts(sce)[hit, , drop = FALSE])
}
for (nm in names(markers)) colData(sce)[[paste0("marker_", nm)]] <- score_module(markers[[nm]])
for (nm in names(programs)) colData(sce)[[paste0("state_score_", nm)]] <- score_module(programs[[nm]])

md <- as.data.table(as.data.frame(colData(sce)), keep.rownames = "cell_id")
mcols <- paste0("marker_", names(markers)); M <- as.matrix(md[, ..mcols])
top <- max.col(M, ties.method = "first"); ord <- t(apply(M, 1, sort, decreasing = TRUE))
md[, coarse_cell_type := names(markers)[top]]
md[, annotation_margin := ord[, 1] - ord[, 2]]
md[, annotation_confidence := pmax(0, pmin(1, annotation_margin / 0.25))]
md[!is.finite(annotation_margin) | annotation_margin < 0.03 | ord[, 1] < 0.02,
   coarse_cell_type := "Unknown_or_ambiguous"]

pcols <- paste0("state_score_", names(programs))
st <- md[coarse_cell_type != "Unknown_or_ambiguous",
  c(list(n_cells = .N, annotation_confidence = mean(annotation_confidence)),
    lapply(.SD, median, na.rm = TRUE)),
  by = .(sample_id, donor, age_class, coarse_cell_type), .SDcols = pcols]
for (v in pcols) {
  mu <- mean(st[[v]], na.rm = TRUE); ss <- sd(st[[v]], na.rm = TRUE)
  st[[paste0("z_", v)]] <- (st[[v]] - mu) / ifelse(ss > 0, ss, 1)
}
zcols <- paste0("z_", pcols); Z <- as.matrix(st[, ..zcols])
Z <- sweep(Z, 1, apply(Z, 1, max), "-") / 0.8; P <- exp(Z); P <- P / rowSums(P)
set(st, j = "phase", value = as.numeric(P %*% c(.12, .28, .45, .64, .80, .94)))
expected <- st[age_class == "young", .(expected_phase = mean(phase),
  young_sd = sd(phase)), by = coarse_cell_type]
st <- merge(st, expected, by = "coarse_cell_type", all.x = TRUE)
st[, deviation := phase - expected_phase]
sigma <- max(0.05, sd(st[age_class == "young"]$deviation, na.rm = TRUE))
rsi <- st[n_cells >= 10 & is.finite(deviation), .(
  RSI = regenerative_synchrony_index(deviation,
    pmin(1, annotation_confidence * (1 - exp(-n_cells / 50))), sigma_reference = sigma),
  phase_dispersion = weighted_var(deviation,
    pmin(1, annotation_confidence * (1 - exp(-n_cells / 50)))),
  effective_cell_types = .N), by = .(sample_id, donor, age_class)]

test <- rsi[, .(n_donors = .N, mean_RSI = mean(RSI, na.rm = TRUE),
  sd_RSI = sd(RSI, na.rm = TRUE), mean_dispersion = mean(phase_dispersion, na.rm = TRUE)),
  by = age_class]
if (all(c("young", "old") %in% rsi$age_class)) {
  rf <- rsi[is.finite(RSI)]; df <- rsi[is.finite(phase_dispersion)]
  wt <- wilcox.test(RSI ~ age_class, data = rf, exact = FALSE)
  wt_disp <- wilcox.test(phase_dispersion ~ age_class, data = df, exact = FALSE)
  auc <- as.numeric(pROC::auc(rf$age_class, -rf$RSI, levels = c("young", "old"),
                              direction = "<"))
  effect <- data.table(contrast = "old_minus_young",
    mean_RSI_difference = mean(rf[age_class == "old"]$RSI) - mean(rf[age_class == "young"]$RSI),
    RSI_Wilcoxon_p = wt$p.value,
    mean_dispersion_difference = mean(df[age_class == "old"]$phase_dispersion) -
      mean(df[age_class == "young"]$phase_dispersion),
    dispersion_Wilcoxon_p = wt_disp$p.value,
    age_classification_AUC_from_negative_RSI = auc,
    n_RSI_old = rf[age_class == "old", .N], n_RSI_young = rf[age_class == "young", .N])
} else effect <- data.table()

qc_summary <- rbindlist(qc_rows); thresholds <- rbindlist(threshold_rows)
all_qc <- rbindlist(all_qc_rows)
fwrite(qc_summary, file.path(outdir, "qc_summary.tsv"), sep = "\t")
fwrite(thresholds, file.path(outdir, "qc_thresholds.tsv"), sep = "\t")
fwrite(all_qc, file.path(outdir, "all_cell_QC_flags.tsv.gz"), sep = "\t")
fwrite(md, file.path(outdir, "cell_metadata.tsv.gz"), sep = "\t")
fwrite(st, file.path(outdir, "sample_celltype_phase.tsv"), sep = "\t")
fwrite(rsi, file.path(outdir, "donor_RSI.tsv"), sep = "\t")
fwrite(test, file.path(outdir, "age_RSI_summary.tsv"), sep = "\t")
fwrite(effect, file.path(outdir, "age_RSI_effect.tsv"), sep = "\t")
writeLines(c(
  "raw counts: public GSE167186 10x H5 files",
  "QC thresholds: donor-specific robust distributions of detected genes, UMIs and mitochondrial fraction",
  "doublets: scDblFinder run independently per donor",
  "annotation: conservative human skeletal-muscle marker modules with ambiguity threshold",
  "inference unit: donor; cells are measurement units only",
  paste("reference sigma:", sigma)
), file.path(outdir, "analysis_manifest.txt"))

p <- ggplot(rsi[is.finite(RSI)], aes(age_class, RSI, colour = age_class)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .08, height = 0) +
  theme_classic(base_size = 9) + guides(colour = "none") + labs(x = NULL, y = "Human muscle RSI")
ggsave(file.path(outdir, "human_age_RSI.pdf"), p, width = 4.2, height = 3.7)

for (v in c("total_counts", "detected_genes", "mitochondrial_percent"))
  set(all_qc, j = v, value = as.numeric(all_qc[[v]]))
qclong <- melt(all_qc, id.vars = "sample_id",
  measure.vars = c("total_counts", "detected_genes", "mitochondrial_percent"),
  variable.name = "metric", value.name = "value")
tl <- rbindlist(list(
  thresholds[, .(sample_id, metric = "total_counts", lower = low_sum, upper = high_sum)],
  thresholds[, .(sample_id, metric = "detected_genes", lower = low_detected, upper = high_detected)],
  thresholds[, .(sample_id, metric = "mitochondrial_percent", lower = 0, upper = high_mt)]
))
pqc <- ggplot(qclong, aes(sample_id, value)) +
  geom_boxplot(outlier.shape = NA, linewidth = .25) +
  geom_point(data = tl, aes(sample_id, lower), colour = "#2166AC", size = .7) +
  geom_point(data = tl, aes(sample_id, upper), colour = "#B2182B", size = .7) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) + theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  labs(x = NULL, y = NULL)
ggsave(file.path(outdir, "QC_threshold_justification.pdf"), pqc, width = 9, height = 8)
cat("GSE167186 human age validation complete:", ncol(sce), "QC singlets from",
    nrow(rsi), "donors.\n")

  },
  "trajectory" = function() {
suppressPackageStartupMessages({
  library(data.table)
  library(slingshot)
  library(FNN)
  library(ggplot2)
})

indir <- "results/complete_dynamics"
outdir <- "results/trajectory_inference"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(indir, "harmonized_sample_celltype_atlas.tsv.gz"))
zcols <- paste0("z_state_score_", c("Activation", "Proliferation", "Commitment",
                                    "Fusion", "Remodeling", "Resolution"))
x[, row_id := .I]
ref <- x[age_class == "young" & injury_model_std %in% c("NTX", "CTX", "BaCl2") &
           complete.cases(x[, ..zcols])]

fit_trajectory <- function(train) {
  Z <- as.matrix(train[, ..zcols]); mu <- colMeans(Z); sc <- apply(Z, 2, sd)
  sc[!is.finite(sc) | sc < 1e-8] <- 1
  Zs <- sweep(sweep(Z, 2, mu, "-"), 2, sc, "/")
  centers <- min(5L, max(3L, floor(nrow(train) / 8L)))
  if (nrow(train) <= centers + 3L) return(NULL)
  set.seed(20260825 + nrow(train))
  cl <- kmeans(Zs, centers = centers, nstart = 30, iter.max = 200)$cluster
  ct <- data.table(cluster = cl, tau = train$target_tau)[, .(tau = mean(tau)), by = cluster]
  start <- as.character(ct$cluster[which.min(ct$tau)])
  end <- as.character(ct$cluster[which.max(ct$tau)])
  ss <- try(slingshot(Zs, clusterLabels = as.character(cl),
                      start.clus = start, end.clus = end), silent = TRUE)
  if (inherits(ss, "try-error")) return(NULL)
  pt <- slingPseudotime(ss); sw <- slingCurveWeights(ss)
  if (is.null(dim(pt))) pt <- matrix(pt, ncol = 1)
  if (is.null(dim(sw))) sw <- matrix(sw, ncol = 1)
  num <- rowSums(pt * sw, na.rm = TRUE); den <- rowSums(sw * is.finite(pt), na.rm = TRUE)
  tau <- num / pmax(den, .Machine$double.eps)
  tau <- (tau - min(tau, na.rm = TRUE)) / diff(range(tau, na.rm = TRUE))
  if (cor(tau, train$target_tau, use = "complete.obs") < 0) tau <- 1 - tau
  list(Z = Zs, tau = tau, mu = mu, sc = sc, clusters = cl,
       start_cluster = start, end_cluster = end)
}

predict_trajectory <- function(fit, test, k = 12L) {
  Z <- as.matrix(test[, ..zcols]); Z <- sweep(sweep(Z, 2, fit$mu, "-"), 2, fit$sc, "/")
  nn <- get.knnx(fit$Z, Z, k = min(k, nrow(fit$Z)))
  w <- 1 / pmax(nn$nn.dist, 1e-4)
  rowSums(matrix(fit$tau[nn$nn.index], nrow = nrow(nn$nn.index)) * w) / rowSums(w)
}

cv_rows <- list()
for (held in unique(ref$dataset_id)) for (ct in unique(ref$cell_type)) {
  tr <- ref[dataset_id != held & cell_type == ct]
  te <- ref[dataset_id == held & cell_type == ct]
  if (nrow(tr) < 18L || nrow(te) < 1L || uniqueN(tr$target_tau) < 4L) next
  fit <- fit_trajectory(tr); if (is.null(fit)) next
  pred <- predict_trajectory(fit, te)
  cv_rows[[paste(held, ct)]] <- data.table(held_out_study = held, cell_type = ct,
    sample_uid = te$sample_uid, target_tau = te$target_tau,
    trajectory_tau = pred, absolute_error = abs(pred - te$target_tau))
}
cv <- rbindlist(cv_rows, fill = TRUE)
cv_sample <- cv[, .(target_tau = unique(target_tau)[1],
  trajectory_tau = median(trajectory_tau, na.rm = TRUE)),
  by = .(held_out_study, sample_uid)]
cv_sample[, absolute_error := abs(trajectory_tau - target_tau)]
cv_metrics <- cv_sample[, .(n_samples = .N, MAE_tau = mean(absolute_error),
  phase_order_correlation = cor(trajectory_tau, target_tau, method = "spearman")),
  by = held_out_study]

# Final reference trajectories and transfer to all ages/injury models.
x[, slingshot_tau := NA_real_]
manifest <- list()
for (ct in unique(ref$cell_type)) {
  tr <- ref[cell_type == ct]
  te <- x[cell_type == ct & complete.cases(x[, ..zcols])]
  if (nrow(tr) < 18L || uniqueN(tr$target_tau) < 4L) next
  fit <- fit_trajectory(tr); if (is.null(fit)) next
  pred <- predict_trajectory(fit, te)
  x$slingshot_tau[match(te$row_id, x$row_id)] <- pred
  manifest[[ct]] <- data.table(cell_type = ct, n_reference_rows = nrow(tr),
    n_clusters = uniqueN(fit$clusters), start_cluster = fit$start_cluster,
    end_cluster = fit$end_cluster,
    reference_spearman = cor(fit$tau, tr$target_tau, method = "spearman"))
}
sample_transfer <- x[is.finite(slingshot_tau), .(
  slingshot_consensus_tau = median(slingshot_tau), n_lineages = .N,
  target_tau = unique(target_tau)[1], time_days = unique(time_days)[1],
  age_class = unique(age_class)[1], injury_model_std = unique(injury_model_std)[1]),
  by = .(dataset_id, sample_id, sample_uid)]

fwrite(cv, file.path(outdir, "leave_study_out_lineage_predictions.tsv.gz"), sep = "\t")
fwrite(cv_sample, file.path(outdir, "leave_study_out_sample_predictions.tsv"), sep = "\t")
fwrite(cv_metrics, file.path(outdir, "leave_study_out_trajectory_metrics.tsv"), sep = "\t")
fwrite(x[, .(row_id, dataset_id, sample_id, sample_uid, cell_type, time_days,
             age_class, injury_model_std, target_tau, slingshot_tau)],
       file.path(outdir, "all_lineage_trajectory_pseudotime.tsv.gz"), sep = "\t")
fwrite(sample_transfer, file.path(outdir, "all_sample_trajectory_clock.tsv"), sep = "\t")
fwrite(rbindlist(manifest), file.path(outdir, "trajectory_fit_manifest.tsv"), sep = "\t")
writeLines(c(
  "trajectory method: slingshot principal curves in six-dimensional sample-lineage state space",
  "orientation: experimental dpi anchors start/end only; curve geometry is learned from expression states",
  "validation: trajectory and nearest-neighbor transfer are refit entirely within each leave-study-out fold",
  "unit: sample-by-lineage pseudobulk state; no cell-level replication"
), file.path(outdir, "analysis_manifest.txt"))

p <- ggplot(cv_sample, aes(target_tau, trajectory_tau, colour = held_out_study)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") + geom_point() +
  theme_classic(base_size = 9) + labs(x = "Anchored biological time", y = "Held-out slingshot time")
ggsave(file.path(outdir, "slingshot_leave_study_out.pdf"), p, width = 5.5, height = 4.2)
cat("Slingshot trajectory inference complete:", nrow(cv_sample), "held-out samples.\n")

  }
)
if(stage=="help") {cat(paste(names(stages),collapse="\n"),"\n")} else {
 if(!stage %in% names(stages)) stop("Unknown stage: ",stage)
 stages[[stage]]()
}
