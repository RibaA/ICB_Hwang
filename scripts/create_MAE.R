#!/usr/bin/env Rscript

# ============================================================
# MultiAssayExperiment creation: GSE136961 (Hwang)
# ============================================================
# Inputs:
#   data/curation/CLIN.csv
#   data/curation/EXPR.csv
#   data/curation/cased_sequenced.csv
#   data/annot/Gencode.v19.annotation.RData
# Output:
#   data/curation/ICB_Hwang.rds
#
# Purpose:
#   Construct an expression-only MultiAssayExperiment using
#   Gencode v19 (hg19) annotations. Preserve GEO clinical
#   metadata and standardized ICB variables in colData.
#   Genes without a matching Gencode v19 symbol are excluded
#   from the annotated MAE but retained in EXPR.csv.
# ============================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else "data/curation"
annot_file <- if (length(args) >= 2) args[2] else "data/annot/Gencode.v19.annotation.RData"
output_file <- if (length(args) >= 3) args[3] else file.path(input_dir, "ICB_Hwang.rds")

load(annot_file) # supplies features_gene

clin <- read.csv(file.path(input_dir, "CLIN.csv"), sep = ";",
                 stringsAsFactors = FALSE, check.names = FALSE)
case <- read.csv(file.path(input_dir, "cased_sequenced.csv"), sep = ";",
                 stringsAsFactors = FALSE)
expr <- read.csv(file.path(input_dir, "EXPR.csv"), sep = ";",
                 row.names = 1, check.names = FALSE)

rownames(clin) <- clin$patient
samples <- intersect(colnames(expr), case$patient[case$expr == 1])
samples <- intersect(samples, rownames(clin))
expr <- as.matrix(expr[, samples, drop = FALSE])
storage.mode(expr) <- "numeric"
clin <- clin[samples, , drop = FALSE]

# Match Hwang gene symbols to v19 gene identifiers (first hit on ambiguity).
feature_idx <- match(rownames(expr), features_gene$gene_name)
keep <- !is.na(feature_idx)
expr <- expr[keep, , drop = FALSE]
feature_rows <- features_gene[feature_idx[keep], , drop = FALSE]

# Match the annotated assay feature identifiers with row ranges.
gene_ids <- as.character(feature_rows$gene_id)
rownames(expr) <- gene_ids
ranges <- GenomicRanges::makeGRangesFromDataFrame(
  feature_rows,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)
names(ranges) <- gene_ids

# Apply the same clinical field renaming convention as ICB_Common.
renamed_cols <- c(
  drug_type = "treatment", primary = "cancer_type",
  t.pfs = "survival_time_pfs", pfs = "event_occurred_pfs",
  t.os = "survival_time_os", os = "event_occurred_os",
  patient = "patientid"
)
for (old in names(renamed_cols)) {
  names(clin)[names(clin) == old] <- renamed_cols[[old]]
}

se <- SummarizedExperiment(
  assays = list(expr = expr),
  rowRanges = ranges,
  colData = S4Vectors::DataFrame(clin)
)
mae <- MultiAssayExperiment(
  experiments = list(expr = se),
  colData = S4Vectors::DataFrame(clin)
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(mae, output_file)
message("Saved MAE: ", output_file)
message("Annotated genes: ", nrow(expr), "; samples: ", ncol(expr))
