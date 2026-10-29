
#!/usr/bin/env Rscript

# ============================================================
# Gene expression harmonization: GSE136961 (Hwang)
# ============================================================
#
# Inputs:
#   data/download/EXP_TPM.tsv
#   data/curation/cased_sequenced.csv
#
# Output:
#   data/curation/EXPR.csv (semicolon-delimited)
#
# Purpose:
#   Harmonize gene expression data for the Hwang cohort
#   following the common ICB expression format.
#
# Processing:
#   1. Read the TPM expression matrix.
#   2. Extract gene symbols from Symbol_ID.
#   3. Resolve duplicate gene symbols by retaining the row
#      with the highest summed TPM across all samples.
#   4. Retain samples with expression data available (expr=1).
#   5. Transform expression values using log2(TPM + 0.001).
#   6. Export a gene-by-sample expression matrix.
#
# Notes:
#   - Input values are TPM, not raw read counts.
#   - Input expression data comes from a targeted immune
#     response gene panel, not whole-transcriptome RNA-seq.
#   - Gene identifiers are harmonized to symbols.
#   - Output rows represent genes and columns represent
#     GEO sample accessions (GSM identifiers).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# 1. Define directories
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_dir  <- if (length(args) >= 1) args[1] else "data/download"
output_dir <- if (length(args) >= 2) args[2] else "data/curation"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read TPM expression matrix
# ------------------------------------------------------------

expr <- data.table::fread(
  file.path(input_dir, "EXP_TPM.tsv"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# Extract gene symbols from Symbol_ID
# Example: ABCF1_11741283 -> ABCF1
gene_symbols <- sub("_.*$", "", expr$Symbol_ID)

# Convert expression columns to numeric matrix
tpm <- as.matrix(
  expr[, setdiff(colnames(expr), "Symbol_ID"), with = FALSE]
)

storage.mode(tpm) <- "numeric"
rownames(tpm) <- gene_symbols

# ------------------------------------------------------------
# 3. Resolve duplicate gene symbols
# ------------------------------------------------------------
# Retain the row with the highest total TPM for each gene
row_totals <- rowSums(tpm, na.rm = TRUE)

# Order by gene symbol and decreasing total expression
ord <- order(rownames(tpm), -row_totals)

tpm <- tpm[ord, , drop = FALSE]

# Keep the first row for each unique gene
tpm <- tpm[
  !duplicated(rownames(tpm)),
  ,
  drop = FALSE
]

# Sort genes alphabetically
tpm <- tpm[
  order(rownames(tpm)),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 4. Select samples with expression data
# ------------------------------------------------------------

case <- read.csv(
  file.path(output_dir, "cased_sequenced.csv"),
  stringsAsFactors = FALSE,
  sep = ";"
)

expr_samples <- case$patient[case$expr == 1]

# Retain only samples with available expression data
tpm <- tpm[
  ,
  colnames(tpm) %in% expr_samples,
  drop = FALSE
]

# ------------------------------------------------------------
# 5. Log2 transformation
# ------------------------------------------------------------

tpm <- log2(tpm + 1)

# ------------------------------------------------------------
# 6. Export harmonized expression matrix
# ------------------------------------------------------------

output_file <- file.path(output_dir, "EXPR.csv")

write.table(
  tpm,
  file = output_file,
  quote = FALSE,
  sep = ";",
  col.names = TRUE,
  row.names = TRUE,
  na = "NA"
)

message("Harmonized expression matrix saved: ", output_file)
message("Genes: ", nrow(tpm))
message("Samples: ", ncol(tpm))
message("Expression harmonization completed successfully.")
