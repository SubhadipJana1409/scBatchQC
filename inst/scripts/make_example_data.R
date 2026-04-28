## ------------------------------------------------------------------
## Script: make_example_data.R
## Purpose: Download, subset, and save real PBMC data from TENxPBMCData
##          for use in inst/extdata/pbmc_small.rds
##
## This script documents the provenance of all data files shipped in
## inst/extdata/, as required by Bioconductor packaging guidelines.
##
## Data source:
##   TENxPBMCData Bioconductor ExperimentHub package
##   - "pbmc3k": ~2,700 PBMCs from a healthy donor (10x Genomics)
##   - "pbmc4k": ~4,340 PBMCs from a healthy donor (10x Genomics)
##   Reference: Zheng GXY et al. (2017). Massively parallel digital
##              transcriptional profiling of single cells.
##              Nature Communications, 8:14049.
##              doi:10.1038/ncomms14049
##
## Processing steps:
##   1. Download PBMC 3k and 4k datasets from TENxPBMCData
##   2. Restrict to shared genes (Ensembl IDs)
##   3. Identify MT genes and top variable genes per batch
##   4. Subset each to 150 cells (random) and ~200 genes
##   5. Assign batch labels ("PBMC_3k", "PBMC_4k")
##   6. Build a clean in-memory SCE with gene symbols as rownames
##   7. Save as a compact SingleCellExperiment to inst/extdata/
##
## To regenerate the data, run from the package root:
##   source("inst/scripts/make_example_data.R")
##
## Requirements:
##   BiocManager::install(c("TENxPBMCData", "scuttle",
##                          "SingleCellExperiment"))
## ------------------------------------------------------------------

suppressPackageStartupMessages({
    library(TENxPBMCData)
    library(SingleCellExperiment)
    library(scuttle)
})

set.seed(2026)

## ---- Step 1: Download real PBMC datasets ----
message("Downloading PBMC 3k...")
pbmc3k <- TENxPBMCData("pbmc3k")

message("Downloading PBMC 4k...")
pbmc4k <- TENxPBMCData("pbmc4k")

## ---- Step 2: Harmonise gene IDs ----
## Both datasets use Ensembl IDs in rownames; keep only shared genes
shared_genes <- intersect(rownames(pbmc3k), rownames(pbmc4k))
message("Shared genes: ", length(shared_genes))

pbmc3k <- pbmc3k[shared_genes, ]
pbmc4k <- pbmc4k[shared_genes, ]

## ---- Step 3: Identify MT genes ----
## Use Symbol_TENx column from pbmc3k to find mitochondrial genes
symbols <- rowData(pbmc3k)$Symbol_TENx
mt_idx <- grep("^MT-", symbols, ignore.case = FALSE)
message("Mitochondrial genes found: ", length(mt_idx))

## ---- Step 4: Select top variable genes + all MT genes ----
## Compute variance PER BATCH separately to avoid cbind rowData conflict
## Realise counts to dense matrices for variance computation
message("Computing gene variances (this may take a moment)...")
vars_3k <- rowVars(as.matrix(counts(pbmc3k)))
vars_4k <- rowVars(as.matrix(counts(pbmc4k)))

## Average variance across batches
avg_vars <- (vars_3k + vars_4k) / 2
names(avg_vars) <- shared_genes

## Pick top 180 most variable non-MT genes, then add all MT genes
non_mt_idx <- setdiff(seq_along(avg_vars), mt_idx)
top_var_order <- order(avg_vars[non_mt_idx], decreasing = TRUE)
n_non_mt <- min(180, length(top_var_order))
keep_genes <- union(
    shared_genes[non_mt_idx[top_var_order[seq_len(n_non_mt)]]],
    shared_genes[mt_idx]
)
message("Genes kept (variable + MT): ", length(keep_genes))

## ---- Step 5: Subset cells ----
## Randomly sample 150 cells from each batch
n_per_batch <- 150

idx_3k <- sample(ncol(pbmc3k), min(n_per_batch, ncol(pbmc3k)))
idx_4k <- sample(ncol(pbmc4k), min(n_per_batch, ncol(pbmc4k)))

sub3k <- pbmc3k[keep_genes, idx_3k]
sub4k <- pbmc4k[keep_genes, idx_4k]

## ---- Step 6: Build clean in-memory SCE ----
## Use gene symbols as rownames for readability
gene_symbols <- rowData(sub3k)$Symbol_TENx

## Realise HDF5-backed counts to in-memory dense matrices
mat3k <- as.matrix(counts(sub3k))
mat4k <- as.matrix(counts(sub4k))

## Combine count matrices directly (no SCE cbind needed)
combined_counts <- cbind(mat3k, mat4k)
rownames(combined_counts) <- gene_symbols
colnames(combined_counts) <- paste0("Cell", seq_len(ncol(combined_counts)))

## Build the final SCE
sce <- SingleCellExperiment(
    assays = list(counts = combined_counts)
)
sce$batch <- rep(c("PBMC_3k", "PBMC_4k"), c(length(idx_3k), length(idx_4k)))

## ---- Step 7: Save ----
outdir <- file.path("inst", "extdata")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile <- file.path(outdir, "pbmc_small.rds")
saveRDS(sce, file = outfile, compress = "xz")

fsize <- file.size(outfile)
message(
    "\n=== Data saved ===\n",
    "File:       ", outfile, "\n",
    "Size:       ", round(fsize / 1024, 1), " KB\n",
    "Dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells\n",
    "Batches:    ", paste(unique(sce$batch), collapse = ", "), "\n",
    "MT genes:   ", sum(grepl("^MT-", rownames(sce))), "\n"
)
