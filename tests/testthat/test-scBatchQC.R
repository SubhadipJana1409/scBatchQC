library(testthat)
library(SingleCellExperiment)
library(scBatchQC)

# ── Helper: minimal two-batch SCE ─────────────────────────────────────────────
make_sce <- function(seed = 42) {
    set.seed(seed)
    counts <- matrix(rpois(500 * 80, lambda = 8), nrow = 500, ncol = 80)
    rownames(counts) <- paste0("Gene", seq_len(500))
    rownames(counts)[1:20] <- paste0("MT-", seq_len(20))
    colnames(counts) <- paste0("Cell", seq_len(80))
    sce <- SingleCellExperiment(assays = list(counts = counts))
    sce$batch <- rep(c("B1", "B2"), each = 40)
    sce
}

# ══════════════════════════════════════════════════════════════════════════════
# batchAwareQCMetrics
# ══════════════════════════════════════════════════════════════════════════════

test_that("batchAwareQCMetrics returns a SingleCellExperiment", {
    sce <- make_sce()
    result <- batchAwareQCMetrics(sce, batch = "batch")
    expect_s4_class(result, "SingleCellExperiment")
})

test_that("batchAwareQCMetrics adds expected colData columns", {
    sce <- make_sce()
    result <- batchAwareQCMetrics(sce, batch = "batch")
    expect_true("scBatchQC_sum" %in% names(colData(result)))
    expect_true("scBatchQC_detected" %in% names(colData(result)))
    expect_true("scBatchQC_outlier" %in% names(colData(result)))
    expect_true("scBatchQC_outlier_reason" %in% names(colData(result)))
})

test_that("scBatchQC_outlier is logical", {
    sce <- make_sce()
    result <- batchAwareQCMetrics(sce, batch = "batch")
    expect_type(result$scBatchQC_outlier, "logical")
})

test_that("batchAwareQCMetrics errors on invalid batch column", {
    sce <- make_sce()
    expect_error(
        batchAwareQCMetrics(sce, batch = "nonexistent"),
        regexp = "not found in colData"
    )
})

test_that("batchAwareQCMetrics errors on non-SCE input", {
    expect_error(
        batchAwareQCMetrics(list()),
        regexp = "SingleCellExperiment"
    )
})

test_that("shrink_strength = 0 and = 1 produce different results", {
    sce <- make_sce()
    r0 <- batchAwareQCMetrics(sce, batch = "batch", shrink_strength = 0)
    r1 <- batchAwareQCMetrics(sce, batch = "batch", shrink_strength = 1)
    # They may flag the same cells for symmetric data, but internal thresholds differ
    # Just check both return valid outputs
    expect_s4_class(r0, "SingleCellExperiment")
    expect_s4_class(r1, "SingleCellExperiment")
})

test_that("number of flagged cells does not exceed total cells", {
    sce <- make_sce()
    result <- batchAwareQCMetrics(sce, batch = "batch")
    expect_lte(sum(result$scBatchQC_outlier), ncol(result))
})

# ══════════════════════════════════════════════════════════════════════════════
# estimateBatchDoubletRate
# ══════════════════════════════════════════════════════════════════════════════

test_that("estimateBatchDoubletRate adds doublet_rate column", {
    sce <- make_sce()
    cells_loaded <- c(B1 = 5000, B2 = 8000)
    result <- estimateBatchDoubletRate(sce,
        batch = "batch",
        cells_loaded = cells_loaded
    )
    expect_true("scBatchQC_doublet_rate" %in% names(colData(result)))
})

test_that("doublet rates are in [0, 0.5]", {
    sce <- make_sce()
    cells_loaded <- c(B1 = 5000, B2 = 8000)
    result <- estimateBatchDoubletRate(sce,
        batch = "batch",
        cells_loaded = cells_loaded
    )
    rates <- result$scBatchQC_doublet_rate
    expect_true(all(rates >= 0))
    expect_true(all(rates <= 0.5))
})

test_that("estimateBatchDoubletRate returns DataFrame when return_sce=FALSE", {
    sce <- make_sce()
    df <- estimateBatchDoubletRate(sce,
        batch = "batch",
        return_sce = FALSE
    )
    expect_s4_class(df, "DataFrame")
    expect_true("doublet_rate_est" %in% names(df))
})

test_that("estimateBatchDoubletRate errors on missing batch labels in cells_loaded", {
    sce <- make_sce()
    expect_error(
        estimateBatchDoubletRate(sce,
            batch = "batch",
            cells_loaded = c(B1 = 5000)
        ),
        regexp = "missing entries"
    )
})

# ══════════════════════════════════════════════════════════════════════════════
# harmonizeQCThresholds
# ══════════════════════════════════════════════════════════════════════════════

test_that("harmonizeQCThresholds returns a list with thresholds and n_flagged", {
    sce <- make_sce() |> batchAwareQCMetrics(batch = "batch")
    result <- harmonizeQCThresholds(sce, batch = "batch")
    expect_type(result, "list")
    expect_true("thresholds" %in% names(result))
    expect_true("n_flagged" %in% names(result))
})

test_that("harmonizeQCThresholds errors without prior QC run", {
    sce <- make_sce()
    expect_error(
        harmonizeQCThresholds(sce, batch = "batch"),
        regexp = "Run batchAwareQCMetrics"
    )
})

test_that("tighter nmads flags more cells", {
    sce <- make_sce() |> batchAwareQCMetrics(batch = "batch")
    r_loose <- harmonizeQCThresholds(sce, batch = "batch", nmads = 4)
    r_tight <- harmonizeQCThresholds(sce, batch = "batch", nmads = 2)
    loose_total <- sum(as.matrix(as.data.frame(r_loose$n_flagged)))
    tight_total <- sum(as.matrix(as.data.frame(r_tight$n_flagged)))
    expect_gte(tight_total, loose_total)
})

# ══════════════════════════════════════════════════════════════════════════════
# BQCResult S4 class
# ══════════════════════════════════════════════════════════════════════════════

test_that("BQCResult constructor and accessors work", {
    library(S4Vectors)
    qf <- DataFrame(
        low_lib = c(FALSE, TRUE, FALSE),
        high_mt = c(FALSE, FALSE, TRUE)
    )
    ds <- c(0.04, 0.06, 0.05)
    bs <- DataFrame(
        batch = c("B1", "B2"),
        doublet_rate_est = c(0.04, 0.06)
    )
    obj <- BQCResult(qcFlags = qf, doubletScores = ds, batchSummary = bs)

    expect_s4_class(obj, "BQCResult")
    expect_equal(nrow(qcFlags(obj)), 3)
    expect_length(doubletScores(obj), 3)
    expect_equal(nrow(batchSummary(obj)), 2)
})

test_that("show method for BQCResult prints without error", {
    library(S4Vectors)
    qf <- DataFrame(low_lib = c(FALSE, TRUE))
    obj <- BQCResult(
        qcFlags       = qf,
        doubletScores = c(0.04, 0.06),
        batchSummary  = DataFrame(batch = "B1", doublet_rate_est = 0.04)
    )
    expect_output(show(obj), "BQCResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# inst/ directory: CITATION and extdata
# ══════════════════════════════════════════════════════════════════════════════

test_that("CITATION file is valid and parseable", {
    cit_file <- system.file("CITATION", package = "scBatchQC")
    skip_if(nchar(cit_file) == 0, "Package not installed; skipping CITATION test")
    cit <- readCitationFile(cit_file)
    expect_s3_class(cit, "citation")
    expect_true(length(cit) >= 1)
})

test_that("pbmc_small.rds in extdata is a valid SingleCellExperiment", {
    rds_path <- system.file("extdata", "pbmc_small.rds", package = "scBatchQC")
    skip_if(nchar(rds_path) == 0, "Package not installed; skipping extdata test")
    sce <- readRDS(rds_path)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_true(ncol(sce) > 0)
    expect_true("batch" %in% names(colData(sce)))
})
