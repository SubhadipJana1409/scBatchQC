# scBatchQC 0.99.2

## Bug fixes

* Improve overall stability.

# scBatchQC 0.99.1

## Bug fixes

* Updated R version dependency from 4.5.0 to 4.6.0 to match Bioconductor 3.23.
* Replaced deprecated `scuttle::perCellQCMetrics()` usage with
  `scrapper::computeRnaQcMetrics()`, removing deprecation warnings from the
  test suite while preserving existing `scBatchQC_*` output column names.
* Updated generated imports and dependency metadata for the new `scrapper`
  QC metric backend.
  
# scBatchQC 0.99.0

## New features

* Initial release submitted to Bioconductor.
* `batchAwareQCMetrics()`: compute per-cell QC metrics with hierarchical
  empirical Bayes batch harmonization.
* `estimateBatchDoubletRate()`: estimate per-batch doublet rates from
  cells-loaded metadata and protocol type.
* `harmonizeQCThresholds()`: inspect and update QC thresholds at arbitrary
  MAD stringency without re-running the full pipeline.
* `plotBatchQC()`: violin plot panel of QC metric distributions per batch
  with harmonized threshold overlays.
* `BQCResult` S4 class for storing batch QC results.
