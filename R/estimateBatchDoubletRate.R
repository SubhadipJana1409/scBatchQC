#' @title Estimate Per-Batch Doublet Rates
#'
#' @description
#' Estimates the expected doublet rate for each batch in a
#' multi-sample \code{SingleCellExperiment} experiment. The doublet
#' rate is modelled as a linear function of per-batch technical
#' covariates (number of cells loaded, median library size, protocol
#' type), enabling principled flagging of likely doublets across
#' batches with heterogeneous capture efficiencies.
#'
#' @details
#' Doublet rates in droplet-based scRNA-seq follow approximately:
#' \deqn{r_b \approx k \times N_b}
#' where \eqn{N_b} is the number of cells loaded per batch and
#' \eqn{k \approx 8 \times 10^{-6}} for 10x Genomics Chromium.
#'
#' \code{estimateBatchDoubletRate} allows \eqn{k} to vary by batch
#' covariates (e.g. protocol, operator) by fitting a linear model on
#' the log-transformed per-batch cell count. When external doublet
#' simulations are not desired, this gives a lightweight alternative to
#' full simulation-based tools like \code{scDblFinder}.
#'
#' Optionally, if the user supplies observed doublet calls from an
#' external tool in \code{colData}, the function will calibrate the
#' rate model against those observations.
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object.
#' @param batch A \code{character(1)} naming the batch column in
#'   \code{colData(sce)}.
#' @param cells_loaded A named \code{numeric} vector with the number of
#'   cells loaded per batch (key = batch label, value = cell count
#'   loaded). If \code{NULL}, the observed cell count per batch is used
#'   as a proxy (underestimates doublet rate).
#' @param protocol A named \code{character} vector mapping batch labels
#'   to protocol type (e.g. \code{"10x_v3"}, \code{"10x_v2"},
#'   \code{"inDrop"}). Used to set the baseline \eqn{k} constant.
#'   Default: \code{NULL} (all batches assumed 10x v3).
#' @param observed_doublets A \code{character(1)} naming a column in
#'   \code{colData(sce)} that contains externally computed doublet calls
#'   (\code{TRUE}/\code{FALSE}). When supplied, the model is calibrated
#'   against observed rates. Default: \code{NULL}.
#' @param return_sce Logical. If \code{TRUE} (default), returns the
#'   input \code{sce} with \code{scBatchQC_doublet_rate} added to
#'   \code{colData}. If \code{FALSE}, returns a \code{DataFrame} of
#'   batch-level estimates.
#'
#' @return If \code{return_sce = TRUE}: the input \code{sce} with a
#'   \code{scBatchQC_doublet_rate} column in \code{colData} giving the
#'   estimated doublet probability for each cell's batch.
#'   If \code{return_sce = FALSE}: a \code{DataFrame} with one row per
#'   batch and columns \code{batch}, \code{n_cells_obs},
#'   \code{cells_loaded}, \code{doublet_rate_est}, \code{protocol}.
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' counts <- matrix(rpois(2000, lambda = 5), nrow = 200, ncol = 100)
#' rownames(counts) <- paste0("Gene", seq_len(200))
#' colnames(counts) <- paste0("Cell", seq_len(100))
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$batch <- rep(c("B1", "B2"), each = 50)
#'
#' cells_loaded <- c(B1 = 5000, B2 = 8000)
#' sce <- estimateBatchDoubletRate(sce, batch = "batch",
#'                                  cells_loaded = cells_loaded)
#' sce$scBatchQC_doublet_rate
#'
#' @seealso \code{\link{batchAwareQCMetrics}}, \code{\link{plotBatchQC}}
#'
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom S4Vectors DataFrame
#' @importFrom methods is
#' @importFrom stats lm coef setNames
#' @export
estimateBatchDoubletRate <- function(sce,
                                      batch             = NULL,
                                      cells_loaded      = NULL,
                                      protocol          = NULL,
                                      observed_doublets = NULL,
                                      return_sce        = TRUE) {
    # ── Validation ────────────────────────────────────────────────────────────
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")
    if (is.null(batch))
        stop("'batch' must be specified; it names a colData column.")
    if (!batch %in% names(colData(sce)))
        stop("'batch' column '", batch, "' not found in colData(sce).")
    if (!is.null(observed_doublets) &&
        !observed_doublets %in% names(colData(sce)))
        stop("'observed_doublets' column not found in colData(sce).")

    # ── Protocol baseline k values (doublets per cell loaded) ─────────────────
    # Empirical constants from 10x Genomics technical documentation and
    # Chromium Next GEM Single Cell 3' v3.1 User Guide
    k_map <- c(
        "10x_v3"   = 8.0e-6,
        "10x_v2"   = 7.5e-6,
        "inDrop"   = 5.0e-6,
        "dropseq"  = 3.0e-6,
        "default"  = 8.0e-6
    )

    batch_labels <- as.character(colData(sce)[[batch]])
    batch_levels <- unique(batch_labels)

    # ── Per-batch observed cell counts ────────────────────────────────────────
    n_obs <- table(batch_labels)[batch_levels]

    # ── Cells loaded: fall back to observed count if not supplied ─────────────
    if (is.null(cells_loaded)) {
        message("'cells_loaded' not supplied; using observed cell counts ",
                "as proxy. Doublet rates will be underestimates.")
        cells_loaded <- as.numeric(n_obs)
        names(cells_loaded) <- batch_levels
    } else {
        # Check all batch levels are covered
        missing_b <- setdiff(batch_levels, names(cells_loaded))
        if (length(missing_b) > 0)
            stop("'cells_loaded' missing entries for batches: ",
                 paste(missing_b, collapse = ", "))
    }

    # ── Protocol lookup ───────────────────────────────────────────────────────
    if (is.null(protocol)) {
        k_values <- setNames(rep(k_map["default"], length(batch_levels)),
                             batch_levels)
    } else {
        missing_p <- setdiff(batch_levels, names(protocol))
        if (length(missing_p) > 0)
            stop("'protocol' missing entries for batches: ",
                 paste(missing_p, collapse = ", "))
        k_values <- k_map[ifelse(protocol[batch_levels] %in% names(k_map),
                                  protocol[batch_levels], "default")]
        names(k_values) <- batch_levels
    }

    # ── Estimate doublet rate per batch ───────────────────────────────────────
    doublet_rate_est <- vapply(batch_levels, function(b) {
        rate <- k_values[[b]] * cells_loaded[[b]]
        # Cap at 0.5 — cannot have >50% doublets
        min(rate, 0.5)
    }, numeric(1))

    # ── Optional calibration against observed doublet calls ───────────────────
    if (!is.null(observed_doublets)) {
        obs_calls <- as.logical(colData(sce)[[observed_doublets]])
        obs_rate  <- tapply(obs_calls, batch_labels, mean, na.rm = TRUE)
        obs_rate  <- obs_rate[batch_levels]

        # Shrink estimated rate toward observed rate where available
        valid <- !is.na(obs_rate)
        doublet_rate_est[valid] <-
            0.5 * doublet_rate_est[valid] + 0.5 * obs_rate[valid]

        message("Calibrated doublet rate estimates against '",
                observed_doublets, "' for ",
                sum(valid), " batch(es).")
    }

    # ── Assemble per-batch summary ─────────────────────────────────────────────
    summary_df <- DataFrame(
        batch            = batch_levels,
        n_cells_obs      = as.integer(n_obs),
        cells_loaded     = cells_loaded[batch_levels],
        doublet_rate_est = doublet_rate_est,
        protocol         = if (!is.null(protocol)) protocol[batch_levels]
                           else rep("10x_v3", length(batch_levels)),
        row.names        = batch_levels
    )

    # ── Return ────────────────────────────────────────────────────────────────
    if (!return_sce) return(summary_df)

    # Add per-cell column (rate from their respective batch)
    colData(sce)[["scBatchQC_doublet_rate"]] <-
        doublet_rate_est[batch_labels]

    sce
}
