#' scBatchQC: Batch-Aware Cell Quality Control for scRNA-seq
#'
#' @description
#' scBatchQC provides a hierarchical empirical Bayes framework for
#' quality control (QC) in multi-sample, multi-batch single-cell
#' RNA-sequencing (scRNA-seq) experiments.
#'
#' Unlike per-sample QC tools such as \code{scuttle::isOutlier}, which
#' apply a single global MAD threshold, \code{scBatchQC} jointly models
#' QC metric distributions (library size, gene count, mitochondrial
#' fraction) and doublet rates across batches. This prevents
#' over-filtering of high-quality batches and under-filtering of
#' low-quality ones — a common but underappreciated problem in
#' multi-batch scRNA-seq workflows.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{batchAwareQCMetrics}}}{Compute per-cell QC
#'     metrics and flag outliers using batch-harmonized MAD thresholds.}
#'   \item{\code{\link{estimateBatchDoubletRate}}}{Model expected
#'     doublet rates per batch as a function of cells loaded and
#'     protocol.}
#'   \item{\code{\link{harmonizeQCThresholds}}}{Inspect and update
#'     harmonized thresholds at arbitrary MAD stringency.}
#'   \item{\code{\link{plotBatchQC}}}{Visualize QC metric distributions
#'     per batch with threshold overlays.}
#' }
#'
#' @section Bioconductor data structures:
#' All functions accept and return
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} objects.
#' Results are stored as additional columns in \code{colData()}.
#'
#' @references
#' Amezquita RA et al. (2020). Orchestrating single-cell analysis with
#' Bioconductor. \emph{Nature Methods}, 17, 137-145.
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rlang .data
#'
#' @docType package
#' @name scBatchQC-package
#' @aliases scBatchQC
"_PACKAGE"
