# scBatchQC

<!-- badges -->
[![Bioc-check](https://github.com/SubhadipJana1409/scBatchQC/actions/workflows/bioc-check.yml/badge.svg)](https://github.com/SubhadipJana1409/scBatchQC/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**scBatchQC** provides a hierarchical empirical Bayes framework
for quality control in multi-sample, multi-batch single-cell
RNA-seq (scRNA-seq) experiments.

Most existing QC tools apply a single global MAD-based threshold
across all cells — ignoring the fact that library size
distributions, mitochondrial fractions, and doublet rates differ
systematically across batches. This leads to **over-filtering of
low-depth batches** and **under-filtering of high-depth ones**.

`scBatchQC` solves this by:

- Estimating per-batch QC metric distributions (median + MAD)
- Shrinking per-batch estimates toward a global empirical Bayes
  prior
- Modeling per-batch doublet rates from cells-loaded metadata
  and protocol type
- Returning a `SingleCellExperiment` with calibrated QC flags
  in `colData`

## The problem scBatchQC solves

```
Standard approach               scBatchQC approach
(scuttle::isOutlier)
────────────────────             ──────────────────
Batch A (fresh, v3)  ──┐        Batch A ← shrunk threshold
Batch B (cryo, v2)   ──┼─ one   Batch B ← shrunk threshold
Batch C (low depth)  ──┘  cut   Batch C ← shrunk threshold

Over-filters Batch C             Calibrated per-batch flags
Under-filters Batch A            Borrow strength across batches
```

## Installation

```r
# Install from Bioconductor (once accepted)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("scBatchQC")

# Or install the development version from GitHub
BiocManager::install("SubhadipJana1409/scBatchQC")
```

## Quick start

```r
library(scBatchQC)
library(SingleCellExperiment)

# Your multi-batch SCE object
# sce$batch must exist in colData

# Step 1: batch-aware QC flagging
sce <- batchAwareQCMetrics(sce, batch = "batch", nmads = 3)

# Step 2: estimate doublet rates
cells_loaded <- c(Batch1 = 8000, Batch2 = 5000)
sce <- estimateBatchDoubletRate(
    sce,
    batch = "batch",
    cells_loaded = cells_loaded
)

# Step 3: visualise
plotBatchQC(sce, batch = "batch")

# Step 4: filter
sce_filtered <- sce[, !sce$scBatchQC_outlier]
```

## Key functions

| Function | Description |
|---|---|
| `batchAwareQCMetrics()` | Per-cell QC with hierarchical MAD thresholds |
| `estimateBatchDoubletRate()` | Per-batch doublet rate modelling |
| `harmonizeQCThresholds()` | Interactive threshold exploration |
| `plotBatchQC()` | QC distribution violin plots per batch |
| `BQCResult` | S4 container for batch QC results |

## How the shrinkage works

For each QC metric and batch, `scBatchQC` estimates:

1. Per-batch median and MAD
2. A global prior pooled across batches (weighted by √n)
3. A shrinkage-adjusted threshold:

```
threshold_b = shrunk_median_b + nmads × shrunk_MAD_b
```

where `shrunk = (1 - s) × per_batch + s × global_prior`
and `s` is the `shrink_strength` parameter (default 0.5).

| `shrink_strength` | Behaviour |
|---|---|
| 0 | Fully per-batch (no pooling) |
| 0.5 | Balanced (default) |
| 1 | Fully pooled (ignores batch) |

## Example output

Tested on real TENxPBMCData (7,040 cells, 2 batches):

```
Outlier cells per batch:
        PBMC_3k PBMC_4k
FALSE      2639    3947
TRUE         61     393

Estimated doublet rates:
PBMC_3k PBMC_4k
  0.040   0.064

Flag rate: 6.4%
```

## Citation

If you use `scBatchQC` in your research, please cite:

> Jana S (2026). scBatchQC: Batch-Aware Cell Quality Control
> for Single-Cell RNA-seq. R package version 0.99.0.
> https://github.com/SubhadipJana1409/scBatchQC

## Contributing

Bug reports and feature requests are welcome at
https://github.com/SubhadipJana1409/scBatchQC/issues

## License

MIT © Subhadip Jana
