#' @title BQCResult: Batch-Aware QC Result Container
#'
#' @description
#' An S4 class to store the output of batch-aware quality control
#' applied to a \code{SingleCellExperiment} object. Slots hold
#' per-cell QC flags, batch-level doublet rate estimates, and
#' harmonized QC thresholds for each batch.
#'
#' @slot qcFlags A \code{DataFrame} with per-cell logical QC flags,
#'   one row per cell and one column per QC metric.
#' @slot doubletScores A \code{numeric} vector of batch-adjusted
#'   doublet probability scores, one per cell.
#' @slot batchSummary A \code{DataFrame} with one row per batch
#'   summarising estimated doublet rate and harmonized thresholds.
#' @slot params A \code{list} storing the parameters used in the
#'   analysis (batch variable name, NMAD multiplier, etc.).
#'
#' @exportClass BQCResult
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
setClass(
    "BQCResult",
    representation(
        qcFlags     = "DataFrame",
        doubletScores = "numeric",
        batchSummary  = "DataFrame",
        params        = "list"
    )
)

#' @title Constructor for BQCResult
#'
#' @description Create a new \code{BQCResult} object.
#'
#' @param qcFlags A \code{DataFrame} of per-cell QC flags.
#' @param doubletScores A numeric vector of doublet scores.
#' @param batchSummary A \code{DataFrame} of per-batch statistics.
#' @param params A list of analysis parameters.
#'
#' @return A \code{BQCResult} object.
#' @export
BQCResult <- function(qcFlags, doubletScores, batchSummary, params = list()) {
    new("BQCResult",
        qcFlags       = qcFlags,
        doubletScores = doubletScores,
        batchSummary  = batchSummary,
        params        = params)
}

# ── Generics ──────────────────────────────────────────────────────────────────

#' @export
setGeneric("qcFlags",
    function(x, ...) standardGeneric("qcFlags"))

#' @export
setGeneric("doubletScores",
    function(x, ...) standardGeneric("doubletScores"))

#' @export
setGeneric("batchSummary",
    function(x, ...) standardGeneric("batchSummary"))

# ── Methods ───────────────────────────────────────────────────────────────────

#' @rdname BQCResult
#' @export
setMethod("show", "BQCResult", function(object) {
    cat("BQCResult\n")
    cat("  Cells          :", nrow(object@qcFlags), "\n")
    cat("  Batches        :", nrow(object@batchSummary), "\n")
    cat("  QC metrics     :", paste(names(object@qcFlags), collapse = ", "), "\n")
    cat("  Doublet range  :",
        round(min(object@doubletScores), 3), "-",
        round(max(object@doubletScores), 3), "\n")
    invisible(object)
})

#' @rdname BQCResult
#' @export
setMethod("qcFlags", "BQCResult", function(x, ...) x@qcFlags)

#' @rdname BQCResult
#' @export
setMethod("doubletScores", "BQCResult", function(x, ...) x@doubletScores)

#' @rdname BQCResult
#' @export
setMethod("batchSummary", "BQCResult", function(x, ...) x@batchSummary)
