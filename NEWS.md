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
