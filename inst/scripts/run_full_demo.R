## ------------------------------------------------------------------
## Script: run_full_demo.R
## Purpose: End-to-end demonstration of all scBatchQC functions
##          using REAL PBMC data from 10x Genomics (Zheng et al. 2017)
##
## Data source:
##   inst/extdata/pbmc_small.rds — a subset of real PBMC 3k + 4k data
##   from TENxPBMCData (Bioconductor ExperimentHub).
##   See inst/scripts/make_example_data.R for full provenance.
##
##   Original reference:
##     Zheng GXY et al. (2017). Massively parallel digital
##     transcriptional profiling of single cells.
##     Nature Communications, 8:14049. doi:10.1038/ncomms14049
##
## This script exercises every exported function in scBatchQC:
##   1. batchAwareQCMetrics()     — batch-harmonised QC flagging
##   2. estimateBatchDoubletRate() — doublet rate modelling
##   3. harmonizeQCThresholds()   — threshold sweep at multiple nmads
##   4. plotBatchQC()             — violin plots with threshold overlays
##   5. BQCResult S4 class        — structured result container
##
## To run from the package root:
##   source("inst/scripts/run_full_demo.R")
##
## Output:
##   - Console report with per-batch QC summaries
##   - scBatchQC_demo_plot.png saved to working directory
## ------------------------------------------------------------------

library(scBatchQC)
library(SingleCellExperiment)
library(S4Vectors)

cat("===== scBatchQC End-to-End Demo (Real PBMC Data) =====\n\n")

# ── 1. Load real PBMC data ────────────────────────────────────

## Try loading from installed package first, fall back to local path
rds_path <- system.file("extdata", "pbmc_small.rds", package = "scBatchQC")
if (nchar(rds_path) == 0) {
    rds_path <- file.path("inst", "extdata", "pbmc_small.rds")
}
stopifnot(
    "pbmc_small.rds not found. Run make_example_data.R first." =
        file.exists(rds_path)
)

sce <- readRDS(rds_path)

cat(sprintf(
    "Loaded real PBMC data: %d genes x %d cells\n",
    nrow(sce), ncol(sce)
))
cat(sprintf(
    "Batches: %s\n",
    paste(unique(sce$batch), collapse = ", ")
))
cat(sprintf(
    "MT genes: %d\n\n",
    sum(grepl("^MT-", rownames(sce)))
))

# ── 2. batchAwareQCMetrics ────────────────────────────────────
## The core function: computes per-cell QC metrics and flags
## outliers using batch-harmonised empirical Bayes thresholds

cat("── Step 1: Batch-aware QC metrics ──\n")

sce <- batchAwareQCMetrics(
    sce,
    batch           = "batch",
    nmads           = 3,
    mt_pattern      = "^MT-",
    shrink_strength = 0.5
)

cat("\nColumns added to colData:\n")
cat(paste(" ", grep("^scBatchQC", names(colData(sce)), value = TRUE),
          collapse = "\n"), "\n")

cat("\nOutlier cells per batch:\n")
print(table(Outlier = sce$scBatchQC_outlier, Batch = sce$batch))

cat("\nMedian library size per batch:\n")
print(tapply(sce$scBatchQC_sum, sce$batch, median))

cat("\nMedian genes detected per batch:\n")
print(tapply(sce$scBatchQC_detected, sce$batch, median))

## Per-batch outlier breakdown
for (b in unique(sce$batch)) {
    mask <- sce$batch == b
    n_out <- sum(sce$scBatchQC_outlier[mask])
    cat(sprintf(
        "\n  %s: %d / %d cells flagged (%.1f%%)\n",
        b, n_out, sum(mask), 100 * n_out / sum(mask)
    ))
}

# ── 3. estimateBatchDoubletRate ───────────────────────────────
## Models expected doublet rate per batch based on cells loaded
## onto the 10x Chromium controller and protocol version

cat("\n── Step 2: Doublet rate estimation ──\n")

## Realistic cells-loaded values for a PBMC experiment
cells_loaded <- c(
    PBMC_3k = 6000,
    PBMC_4k = 8000
)

protocol_map <- c(
    PBMC_3k = "10x_v2",
    PBMC_4k = "10x_v2"
)

sce <- estimateBatchDoubletRate(
    sce,
    batch        = "batch",
    cells_loaded = cells_loaded,
    protocol     = protocol_map
)

cat("\nEstimated doublet rate per batch:\n")
print(tapply(sce$scBatchQC_doublet_rate, sce$batch, unique))

## Also get the DataFrame summary
doublet_summary <- estimateBatchDoubletRate(
    sce,
    batch        = "batch",
    cells_loaded = cells_loaded,
    protocol     = protocol_map,
    return_sce   = FALSE
)
cat("\nDoublet rate summary table:\n")
print(doublet_summary)

# ── 4. harmonizeQCThresholds ──────────────────────────────────
## Sweep different nmads values to see how stringency affects
## the number of flagged cells — helps users pick the right cutoff

cat("\n── Step 3: Threshold exploration ──\n")

cat("\nThreshold sweep across nmads values:\n")
for (n in c(2.0, 2.5, 3.0, 3.5, 4.0)) {
    r <- harmonizeQCThresholds(sce, batch = "batch", nmads = n)
    total <- sum(as.matrix(as.data.frame(r$n_flagged)), na.rm = TRUE)
    cat(sprintf("  nmads = %.1f -> %d cells flagged total\n", n, total))
}

result <- harmonizeQCThresholds(sce, batch = "batch", nmads = 3)
cat("\nPer-batch thresholds at nmads=3:\n")
print(result$thresholds)
cat("\nFlagged cells per batch per metric (nmads=3):\n")
print(result$n_flagged)

# ── 5. plotBatchQC ────────────────────────────────────────────
## Generates violin plots of QC metric distributions per batch
## with harmonised threshold lines and outlier highlights

cat("\n── Step 4: Plotting QC distributions ──\n")

p <- plotBatchQC(sce, batch = "batch")

output_path <- file.path(getwd(), "scBatchQC_demo_plot.png")
ggplot2::ggsave(output_path, plot = p, width = 8, height = 10, dpi = 150)
cat(sprintf("Plot saved to: %s\n", output_path))

# ── 6. BQCResult S4 class ────────────────────────────────────
## Structured container for programmatic access to QC results

cat("\n── Step 5: BQCResult S4 class ──\n")

qc_flags <- DataFrame(
    low_lib_size = sce$scBatchQC_sum < quantile(sce$scBatchQC_sum, 0.05),
    high_mt      = !is.na(sce$scBatchQC_outlier_reason) &
        grepl("MT", sce$scBatchQC_outlier_reason)
)

result_obj <- BQCResult(
    qcFlags       = qc_flags,
    doubletScores = sce$scBatchQC_doublet_rate,
    batchSummary  = doublet_summary
)

show(result_obj)
cat("\nqcFlags accessor (first 6 cells):\n")
print(head(qcFlags(result_obj)))
cat("\ndoubletScores accessor (first 6 cells):\n")
print(head(doubletScores(result_obj)))
cat("\nbatchSummary accessor:\n")
print(batchSummary(result_obj))

# ── 7. Filter and final report ────────────────────────────────

cat("\n── Step 6: Filter and final report ──\n")

sce_filtered <- sce[, !sce$scBatchQC_outlier]

cat(sprintf("Before filtering : %d cells\n", ncol(sce)))
cat(sprintf("After filtering  : %d cells\n", ncol(sce_filtered)))
cat(sprintf("Cells removed    : %d (%.1f%%)\n",
            ncol(sce) - ncol(sce_filtered),
            100 * (ncol(sce) - ncol(sce_filtered)) / ncol(sce)))

cat("\nRetained cells per batch:\n")
print(table(sce_filtered$batch))

# ── 8. Summary ────────────────────────────────────────────────

cat("\n========================================\n")
cat("  ALL STEPS COMPLETED SUCCESSFULLY\n")
cat("  scBatchQC is working correctly\n")
cat("  Data: Real PBMC 3k + 4k\n")
cat("  (Zheng et al. 2017, Nat Comm)\n")
cat("========================================\n")
cat(sprintf("Plot: %s\n", output_path))
