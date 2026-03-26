#' @title Harmonize QC Thresholds Across Batches
#'
#' @description
#' Given a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' that has already been processed by \code{\link{batchAwareQCMetrics}},
#' \code{harmonizeQCThresholds} returns the per-batch QC threshold
#' table and optionally updates \code{colData} with revised flags at a
#' user-specified stringency.
#'
#' This is useful for interactive threshold exploration or for
#' downstream reporting: instead of re-running the full QC pipeline,
#' the user can sweep \code{nmads} and inspect how the number of
#' flagged cells changes per batch.
#'
#' @param sce A \code{SingleCellExperiment} processed by
#'   \code{\link{batchAwareQCMetrics}}.
#' @param batch A \code{character(1)} naming the batch column.
#' @param nmads A \code{numeric(1)} MAD multiplier. Default: \code{3}.
#' @param shrink_strength A \code{numeric(1)} in \code{[0, 1]}.
#'   Default: \code{0.5}.
#' @param update_sce Logical. If \code{TRUE}, rewrites the
#'   \code{scBatchQC_outlier} and \code{scBatchQC_outlier_reason}
#'   columns in \code{colData(sce)} using the new thresholds.
#'   Default: \code{FALSE}.
#'
#' @return A \code{list} with components:
#'   \describe{
#'     \item{\code{thresholds}}{A named list of per-metric
#'       \code{data.frame}s (rows = batches, columns = \code{lower}
#'       and \code{upper}).}
#'     \item{\code{n_flagged}}{A \code{DataFrame} with one row per
#'       batch and one column per metric showing the number of cells
#'       that would be flagged at these thresholds.}
#'     \item{\code{sce}}{The (possibly updated) \code{sce}, returned
#'       only when \code{update_sce = TRUE}.}
#'   }
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' counts <- matrix(rpois(2000, lambda = 5), nrow = 200, ncol = 100)
#' rownames(counts) <- paste0("Gene", seq_len(200))
#' rownames(counts)[1:10] <- paste0("MT-", seq_len(10))
#' colnames(counts) <- paste0("Cell", seq_len(100))
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$batch <- rep(c("B1", "B2"), each = 50)
#' sce <- batchAwareQCMetrics(sce, batch = "batch")
#'
#' # Explore with 2.5 MADs instead of default 3
#' result <- harmonizeQCThresholds(sce, batch = "batch", nmads = 2.5)
#' result$n_flagged
#'
#' @seealso \code{\link{batchAwareQCMetrics}},
#'   \code{\link{estimateBatchDoubletRate}}, \code{\link{plotBatchQC}}
#'
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @export
harmonizeQCThresholds <- function(sce,
                                   batch           = NULL,
                                   nmads           = 3,
                                   shrink_strength = 0.5,
                                   update_sce      = FALSE) {
    # ── Validation ─────────────────────────────────────────────────────────────
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")

    qc_cols <- grep("^scBatchQC_(sum|detected|subsets_MT_percent)$",
                    names(colData(sce)), value = TRUE)
    if (length(qc_cols) == 0)
        stop("No scBatchQC QC metric columns found in colData(sce). ",
             "Run batchAwareQCMetrics() first.")

    if (!is.null(batch) && !batch %in% names(colData(sce)))
        stop("'batch' column '", batch, "' not found in colData(sce).")

    # ── Reconstruct qc_df from existing colData columns ───────────────────────
    metrics <- sub("^scBatchQC_", "", qc_cols)
    qc_df <- as.data.frame(colData(sce)[qc_cols])
    names(qc_df) <- metrics

    batch_labels <- if (!is.null(batch))
        as.character(colData(sce)[[batch]])
    else
        rep("all", ncol(sce))
    batch_levels <- unique(batch_labels)

    # ── Recompute thresholds ───────────────────────────────────────────────────
    thresholds <- .computeHarmonizedThresholds(
        qc_df        = qc_df,
        metrics      = metrics,
        batch_labels = batch_labels,
        batch_levels = batch_levels,
        nmads        = nmads,
        shrink_strength = shrink_strength
    )

    # ── Count flagged cells per batch per metric ───────────────────────────────
    n_flagged_list <- lapply(metrics, function(metric) {
        thresh <- thresholds[[metric]]
        vals   <- qc_df[[metric]]
        upper  <- thresh[batch_labels, "upper"]
        lower  <- thresh[batch_labels, "lower"]
        flagged <- (vals > upper) |
                   (vals < lower & metric != "subsets_MT_percent")
        tapply(flagged, batch_labels, sum)[batch_levels]
    })
    n_flagged_df <- DataFrame(
        do.call(cbind, n_flagged_list),
        row.names = batch_levels
    )
    names(n_flagged_df) <- metrics

    out <- list(thresholds = thresholds, n_flagged = n_flagged_df)

    # ── Optionally rewrite sce ─────────────────────────────────────────────────
    if (update_sce) {
        outlier_mat <- mapply(function(metric) {
            thresh <- thresholds[[metric]]
            vals   <- qc_df[[metric]]
            upper  <- thresh[batch_labels, "upper"]
            lower  <- thresh[batch_labels, "lower"]
            (vals > upper) | (vals < lower & metric != "subsets_MT_percent")
        }, metrics)

        if (is.null(dim(outlier_mat)))
            outlier_mat <- matrix(outlier_mat, ncol = 1,
                                  dimnames = list(NULL, metrics))

        outlier_any    <- rowSums(outlier_mat, na.rm = TRUE) > 0
        outlier_reason <- apply(outlier_mat, 1, function(x)
            paste(metrics[x], collapse = "; "))
        outlier_reason[outlier_reason == ""] <- NA_character_

        colData(sce)[["scBatchQC_outlier"]]        <- outlier_any
        colData(sce)[["scBatchQC_outlier_reason"]] <- outlier_reason

        out$sce <- sce
    }

    out
}
