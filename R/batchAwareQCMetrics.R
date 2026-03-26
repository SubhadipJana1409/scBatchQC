#' @title Batch-Aware QC Metric Computation
#'
#' @description
#' Computes per-cell quality control metrics and identifies outlier cells
#' using a hierarchical empirical Bayes approach. Unlike
#' \code{scuttle::isOutlier},
#' which applies a single global MAD threshold, \code{batchAwareQCMetrics}
#' estimates batch-specific MAD scales and shrinks them toward a global
#' prior, preventing over-filtering in high-quality batches and
#' under-filtering in low-quality ones.
#'
#' @details
#' For each QC metric \eqn{m} and batch \eqn{b}, the function estimates:
#' \enumerate{
#'   \item Per-batch median \eqn{\mu_b} and MAD \eqn{\sigma_b}.
#'   \item A shrinkage weight \eqn{w_b} based on batch cell count.
#'   \item A global prior \eqn{\mu_0, \sigma_0} pooled across batches.
#'   \item A harmonized threshold
#'         \eqn{\tau_b = \mu_b^* + nmads \times \sigma_b^*}
#'         where \eqn{\mu_b^*} and \eqn{\sigma_b^*} are the
#'         shrinkage estimates.
#' }
#' A cell is flagged as an outlier if any QC metric exceeds its
#' batch-specific harmonized threshold.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object. Must have raw counts in \code{assay(sce, "counts")}.
#' @param batch A \code{character(1)} naming a column in \code{colData(sce)}
#'   that identifies batch membership. If \code{NULL}, falls back to
#'   standard per-dataset QC (equivalent to \code{scuttle::isOutlier}).
#' @param metrics A \code{character} vector of QC metrics to evaluate.
#'   Supported: \code{"sum"} (library size), \code{"detected"}
#'   (genes detected), \code{"subsets_MT_percent"} (mitochondrial
#'   percentage). Default: all three.
#' @param nmads A \code{numeric(1)} number of MADs to use as the
#'   outlier threshold. Default: \code{3}.
#' @param mt_pattern A \code{character(1)} regex passed to
#'   \code{scuttle::perCellQCMetrics} to identify mitochondrial genes.
#'   Default: \code{"^MT-"}.
#' @param shrink_strength A \code{numeric(1)} in \code{[0, 1]} controlling
#'   how much per-batch estimates are shrunk toward the global prior.
#'   \code{0} = no shrinkage (pure per-batch); \code{1} = full pooling.
#'   Default: \code{0.5} (empirical Bayes midpoint).
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   controlling parallelisation. Default: \code{SerialParam()}.
#'
#' @return The input \code{sce} with the following additions to
#'   \code{colData}:
#'   \itemize{
#'     \item \code{scBatchQC_sum}: library size (total UMI count).
#'     \item \code{scBatchQC_detected}: number of detected genes.
#'     \item \code{scBatchQC_subsets_MT_percent}: mitochondrial fraction.
#'     \item \code{scBatchQC_outlier}: logical flag; \code{TRUE} if the cell
#'       fails any QC threshold.
#'     \item \code{scBatchQC_outlier_reason}: character string naming which
#'       metric(s) caused the flag.
#'   }
#'
#' @references
#' Lun ATL et al. (2016). A step-by-step workflow for low-level analysis
#' of single-cell RNA sequencing data with Bioconductor.
#' \emph{F1000Research}, 5, 2122.
#'
#' @seealso \code{\link{estimateBatchDoubletRate}},
#'   \code{\link{harmonizeQCThresholds}}, \code{\link{plotBatchQC}}
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' # Simulate a minimal SCE with two batches
#' set.seed(42)
#' counts <- matrix(rpois(2000, lambda = 5), nrow = 200, ncol = 100)
#' rownames(counts) <- paste0("Gene", seq_len(200))
#' rownames(counts)[1:10] <- paste0("MT-", seq_len(10))
#' colnames(counts) <- paste0("Cell", seq_len(100))
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$batch <- rep(c("B1", "B2"), each = 50)
#'
#' sce <- batchAwareQCMetrics(sce, batch = "batch")
#' table(sce$scBatchQC_outlier, sce$batch)
#'
#' @importFrom scuttle perCellQCMetrics
#' @importFrom SummarizedExperiment colData assay "colData<-"
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom stats median mad
#' @export
batchAwareQCMetrics <- function(sce,
                                batch = NULL,
                                metrics = c(
                                    "sum", "detected",
                                    "subsets_MT_percent"
                                ),
                                nmads = 3,
                                mt_pattern = "^MT-",
                                shrink_strength = 0.5,
                                BPPARAM = SerialParam()) {
    # ── Input validation ────────────────────────────────────
    if (!is(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment object.")
    }
    if (!is.null(batch) && !batch %in% names(colData(sce))) {
        stop("'batch' column '", batch, "' not found in colData(sce).")
    }
    if (!is.numeric(nmads) || nmads <= 0) {
        stop("'nmads' must be a positive number.")
    }
    if (!is.numeric(shrink_strength) ||
        shrink_strength < 0 || shrink_strength > 1) {
        stop("'shrink_strength' must be in [0, 1].")
    }

    # ── Compute per-cell QC metrics ─────────────────────────
    mt_idx <- grep(mt_pattern, rownames(sce))
    subsets <- if (length(mt_idx) > 0) {
        list(MT = mt_idx)
    } else {
        list()
    }

    qc_df <- perCellQCMetrics(sce, subsets = subsets, BPPARAM = BPPARAM)

    # Normalise metric names (handle missing MT when no MT genes found)
    available_metrics <- intersect(metrics, names(qc_df))
    if (length(available_metrics) == 0) {
        stop("None of the requested metrics found in perCellQCMetrics output.")
    }

    # ── Determine batch labels ──────────────────────────────
    if (is.null(batch)) {
        batch_labels <- rep("all", ncol(sce))
        message("No 'batch' supplied; using single-group QC.")
    } else {
        batch_labels <- as.character(colData(sce)[[batch]])
    }
    batch_levels <- unique(batch_labels)

    # ── Estimate per-batch + global prior ────────────────────────────────────
    thresholds <- .computeHarmonizedThresholds(
        qc_df = qc_df,
        metrics = available_metrics,
        batch_labels = batch_labels,
        batch_levels = batch_levels,
        nmads = nmads,
        shrink_strength = shrink_strength
    )

    # ── Flag outlier cells ──────────────────────────────────
    outlier_mat <- mapply(function(metric) {
        thresh_upper <- thresholds[[metric]][batch_labels, "upper"]
        thresh_lower <- thresholds[[metric]][batch_labels, "lower"]
        vals <- qc_df[[metric]]
        (vals > thresh_upper) |
            (vals < thresh_lower & metric != "subsets_MT_percent")
    }, available_metrics)

    if (is.null(dim(outlier_mat))) {
        outlier_mat <- matrix(outlier_mat,
            ncol = 1,
            dimnames = list(NULL, available_metrics)
        )
    }

    outlier_any <- rowSums(outlier_mat, na.rm = TRUE) > 0
    outlier_reason <- apply(outlier_mat, 1, function(x) {
        paste(available_metrics[x], collapse = "; ")
    })
    outlier_reason[outlier_reason == ""] <- NA_character_

    # ── Write back to colData ───────────────────────────────
    for (m in available_metrics) {
        colData(sce)[[paste0("scBatchQC_", m)]] <- qc_df[[m]]
    }
    colData(sce)[["scBatchQC_outlier"]] <- outlier_any
    colData(sce)[["scBatchQC_outlier_reason"]] <- outlier_reason

    sce
}


# ── Internal helpers ────────────────────────────────────────

#' Compute harmonized (shrinkage-based) QC thresholds
#'
#' @keywords internal
#' @noRd
.computeHarmonizedThresholds <- function(qc_df, metrics, batch_labels,
                                         batch_levels, nmads,
                                         shrink_strength) {
    thresholds <- lapply(metrics, function(metric) {
        vals <- qc_df[[metric]]

        # Per-batch location and scale estimates
        batch_stats <- lapply(batch_levels, function(b) {
            idx <- batch_labels == b
            v <- vals[idx]
            v <- v[is.finite(v)]
            list(
                n      = sum(idx),
                median = median(v),
                mad    = max(mad(v), 1e-6) # guard against zero MAD
            )
        })
        names(batch_stats) <- batch_levels

        # Global prior (weighted by sqrt(n) for robustness)
        weights <- vapply(batch_stats, function(s) {
            sqrt(s$n)
        }, numeric(1))
        weights <- weights / sum(weights)

        global_median <- sum(
            weights * vapply(batch_stats, `[[`, numeric(1), "median")
        )
        global_mad <- sum(
            weights * vapply(batch_stats, `[[`, numeric(1), "mad")
        )

        # Shrinkage: batch estimate → global prior
        out <- lapply(batch_levels, function(b) {
            s <- batch_stats[[b]]
            shrunken_median <- (1 - shrink_strength) * s$median +
                shrink_strength * global_median
            shrunken_mad <- (1 - shrink_strength) * s$mad +
                shrink_strength * global_mad
            c(
                upper = shrunken_median + nmads * shrunken_mad,
                lower = shrunken_median - nmads * shrunken_mad
            )
        })
        as.data.frame(do.call(rbind, out), row.names = batch_levels)
    })
    names(thresholds) <- metrics
    thresholds
}
