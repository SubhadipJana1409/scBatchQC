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
#' @importFrom methods new is
setClass(
    "BQCResult",
    representation(
        qcFlags       = "DataFrame",
        doubletScores = "numeric",
        batchSummary  = "DataFrame",
        params        = "list"
    )
)

#' @title Constructor for BQCResult
#' @description Create a new \code{BQCResult} object.
#' @param qcFlags A \code{DataFrame} of per-cell QC flags.
#' @param doubletScores A numeric vector of doublet scores.
#' @param batchSummary A \code{DataFrame} of per-batch statistics.
#' @param params A list of analysis parameters.
#' @return A \code{BQCResult} object.
#' @examples
#' library(S4Vectors)
#' qf <- DataFrame(low_lib = c(FALSE, TRUE, FALSE))
#' obj <- BQCResult(
#'     qcFlags       = qf,
#'     doubletScores = c(0.04, 0.06, 0.05),
#'     batchSummary  = DataFrame(batch = "B1", doublet_rate_est = 0.04)
#' )
#' obj
#' @export
BQCResult <- function(qcFlags, doubletScores, batchSummary, params = list()) {
    new("BQCResult",
        qcFlags       = qcFlags,
        doubletScores = doubletScores,
        batchSummary  = batchSummary,
        params        = params)
}

# ── Generics ────────────────────────────────────────────────

#' @title Accessor for QC flags in a BQCResult
#' @description Returns the per-cell QC flag \code{DataFrame}.
#' @param x A \code{BQCResult} object.
#' @param ... Additional arguments (not used).
#' @return A \code{DataFrame} of per-cell logical QC flags.
#' @examples
#' library(S4Vectors)
#' obj <- BQCResult(
#'     qcFlags       = DataFrame(low_lib = c(FALSE, TRUE)),
#'     doubletScores = c(0.04, 0.06),
#'     batchSummary  = DataFrame(batch = "B1", doublet_rate_est = 0.04)
#' )
#' qcFlags(obj)
#' @export
setGeneric("qcFlags", function(x, ...) standardGeneric("qcFlags"))

#' @title Accessor for doublet scores in a BQCResult
#' @description Returns the per-cell doublet score vector.
#' @param x A \code{BQCResult} object.
#' @param ... Additional arguments (not used).
#' @return A \code{numeric} vector of doublet scores.
#' @examples
#' library(S4Vectors)
#' obj <- BQCResult(
#'     qcFlags       = DataFrame(low_lib = c(FALSE, TRUE)),
#'     doubletScores = c(0.04, 0.06),
#'     batchSummary  = DataFrame(batch = "B1", doublet_rate_est = 0.04)
#' )
#' doubletScores(obj)
#' @export
setGeneric("doubletScores", function(x, ...) standardGeneric("doubletScores"))

#' @title Accessor for batch summary in a BQCResult
#' @description Returns the per-batch summary \code{DataFrame}.
#' @param x A \code{BQCResult} object.
#' @param ... Additional arguments (not used).
#' @return A \code{DataFrame} of per-batch statistics.
#' @examples
#' library(S4Vectors)
#' obj <- BQCResult(
#'     qcFlags       = DataFrame(low_lib = c(FALSE, TRUE)),
#'     doubletScores = c(0.04, 0.06),
#'     batchSummary  = DataFrame(batch = "B1", doublet_rate_est = 0.04)
#' )
#' batchSummary(obj)
#' @export
setGeneric("batchSummary", function(x, ...) standardGeneric("batchSummary"))

# ── Methods ─────────────────────────────────────────────────

#' @title Show method for BQCResult
#' @description Prints a compact summary of a \code{BQCResult} object.
#' @param object A \code{BQCResult} object.
#' @return Invisibly returns \code{object}.
#' @importFrom methods show
#' @export
setMethod("show", "BQCResult", function(object) {
    cat("BQCResult\n")
    cat("  Cells          :", nrow(object@qcFlags), "\n")
    cat("  Batches        :", nrow(object@batchSummary), "\n")
    cat("  QC metrics     :",
        paste(names(object@qcFlags), collapse = ", "), "\n")
    cat("  Doublet range  :",
        round(min(object@doubletScores), 3), "-",
        round(max(object@doubletScores), 3), "\n")
    invisible(object)
})

#' @rdname qcFlags
#' @export
setMethod("qcFlags", "BQCResult", function(x, ...) x@qcFlags)

#' @rdname doubletScores
#' @export
setMethod("doubletScores", "BQCResult", function(x, ...) x@doubletScores)

#' @rdname batchSummary
#' @export
setMethod("batchSummary", "BQCResult", function(x, ...) x@batchSummary)
