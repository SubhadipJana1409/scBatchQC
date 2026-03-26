#' @title Visualize QC Metric Distributions Across Batches
#'
#' @description
#' Produces a panel of violin plots showing per-batch distributions of
#' QC metrics, with harmonized threshold lines overlaid. Useful for
#' inspecting whether \code{\link{batchAwareQCMetrics}} thresholds are
#' sensible and for comparing batch quality visually.
#'
#' @param sce A \code{SingleCellExperiment} processed by
#'   \code{\link{batchAwareQCMetrics}}.
#' @param batch A \code{character(1)} naming the batch column in
#'   \code{colData(sce)}.
#' @param metrics A \code{character} vector of QC metric names to plot
#'   (without the \code{"scBatchQC_"} prefix). Default: all available.
#' @param show_thresholds Logical. If \code{TRUE}, overlays the
#'   batch-specific harmonized upper thresholds as dashed horizontal
#'   lines. Default: \code{TRUE}.
#' @param nmads Passed to \code{\link{harmonizeQCThresholds}} to
#'   recompute thresholds for display. Default: \code{3}.
#' @param colour_by A \code{character(1)} naming a \code{colData}
#'   column to colour cells by (e.g. \code{"scBatchQC_outlier"}).
#'   Default: \code{"scBatchQC_outlier"}.
#' @param point_size Numeric. Jitter point size. Default: \code{0.4}.
#' @param point_alpha Numeric. Jitter point alpha. Default: \code{0.4}.
#'
#' @return A \code{ggplot2} object. Can be modified with standard
#'   \code{ggplot2} functions or saved with \code{ggsave()}.
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
#' plotBatchQC(sce, batch = "batch")
#'
#' @seealso \code{\link{batchAwareQCMetrics}},
#'   \code{\link{harmonizeQCThresholds}}
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter facet_wrap
#'   theme_bw labs geom_hline scale_color_manual theme element_text
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @export
plotBatchQC <- function(sce,
                         batch           = NULL,
                         metrics         = NULL,
                         show_thresholds = TRUE,
                         nmads           = 3,
                         colour_by       = "scBatchQC_outlier",
                         point_size      = 0.4,
                         point_alpha     = 0.4) {
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")

    qc_cols <- grep("^scBatchQC_(sum|detected|subsets_MT_percent)$",
                    names(colData(sce)), value = TRUE)
    if (length(qc_cols) == 0)
        stop("No QC columns found. Run batchAwareQCMetrics() first.")

    available_metrics <- sub("^scBatchQC_", "", qc_cols)
    if (is.null(metrics)) {
        metrics <- available_metrics
    } else {
        metrics <- intersect(metrics, available_metrics)
    }

    cd <- as.data.frame(colData(sce))
    if (!is.null(batch) && batch %in% names(cd)) {
        cd$.__batch__. <- cd[[batch]]
    } else {
        cd$.__batch__. <- "all"
    }

    colour_col <- if (!is.null(colour_by) && colour_by %in% names(cd))
        colour_by else NULL

    # Melt into long format (no reshape2 dependency)
    plot_df <- do.call(rbind, lapply(metrics, function(m) {
        col_name <- paste0("scBatchQC_", m)
        data.frame(
            batch  = cd$.__batch__.,
            metric = m,
            value  = cd[[col_name]],
            colour = if (!is.null(colour_col)) as.character(cd[[colour_col]])
                     else "all",
            stringsAsFactors = FALSE
        )
    }))

    # Threshold lines per metric per batch
    thresh_df <- NULL
    if (show_thresholds) {
        result <- harmonizeQCThresholds(sce, batch = batch, nmads = nmads)
        thresh_df <- do.call(rbind, lapply(metrics, function(m) {
            t <- result$thresholds[[m]]
            data.frame(
                batch  = rownames(t),
                metric = m,
                upper  = t$upper,
                stringsAsFactors = FALSE
            )
        }))
    }

    p <- ggplot(plot_df,
                aes(x = .data[["batch"]], y = .data[["value"]],
                    colour = .data[["colour"]])) +
        geom_violin(fill = NA, colour = "grey60", linewidth = 0.4) +
        geom_jitter(width = 0.15, size = point_size, alpha = point_alpha) +
        facet_wrap(~ metric, scales = "free_y", ncol = 1,
                   strip.position = "left") +
        theme_bw(base_size = 11) +
        theme(
            strip.placement  = "outside",
            strip.background = element_blank(),
            axis.text.x      = element_text(angle = 30, hjust = 1)
        ) +
        labs(x = "Batch", y = NULL,
             colour = colour_col,
             title  = "Batch-aware QC metric distributions")

    if (!is.null(colour_col) &&
        all(c("TRUE", "FALSE") %in% unique(plot_df$colour))) {
        p <- p + scale_color_manual(
            values = c("TRUE" = "#E24B4A", "FALSE" = "#888780"),
            labels = c("TRUE" = "Outlier", "FALSE" = "Pass"))
    }

    if (show_thresholds && !is.null(thresh_df)) {
        # Overlay threshold lines (geom_hline per facet via geom_segment workaround)
        p <- p + geom_hline(
            data     = thresh_df,
            aes(yintercept = .data[["upper"]]),
            linetype = "dashed",
            colour   = "#E24B4A",
            linewidth = 0.5,
            alpha    = 0.7,
            inherit.aes = FALSE
        )
    }

    p
}
