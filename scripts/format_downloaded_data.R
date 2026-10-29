
#!/usr/bin/env Rscript

# ============================================================
# GEO data preprocessing: GSE136961 (Hwang)
# ============================================================
#
# Inputs:
#   data/download/GSE136961_series_matrix.txt.gz
#   data/download/GSE136961_TPM.tsv.gz
#
# Outputs:
#   data/download/CLIN.txt
#   data/download/EXP_TPM.tsv
#
# Purpose:
#   Extract GEO clinical metadata and harmonize expression
#   sample identifiers with GEO accessions.
#
# Notes:
#   - Original GEO files remain unchanged.
#   - GEO metadata is decompressed for GEOquery.
#   - Expression data are read directly from the gzip file.
#   - Expression sample titles are mapped to GEO accessions.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
})

# ------------------------------------------------------------
# 1. Define input directory
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

work_dir <- if (length(args) >= 1) args[1] else "data/download"

# ------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------

clin_gz <- file.path(work_dir, "GSE136961_series_matrix.txt.gz")
expr_gz <- file.path(work_dir, "GSE136961_TPM.tsv.gz")

clin_txt <- file.path(work_dir, "GSE136961_series_matrix.txt")

clin_output <- file.path(work_dir, "CLIN.txt")
expr_output <- file.path(work_dir, "EXP_TPM.tsv")

# ------------------------------------------------------------
# 3. Extract GEO clinical metadata
# ------------------------------------------------------------

message("Reading GEO sample metadata...")

# Decompress GEO metadata if needed
if (!file.exists(clin_txt)) {
  con_in <- gzfile(clin_gz, "rb")
  con_out <- file(clin_txt, "wb")

  tryCatch({
    repeat {
      chunk <- readBin(con_in, "raw", n = 1024 * 1024)
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }
  }, finally = {
    close(con_in)
    close(con_out)
  })
}

geo <- GEOquery::getGEO(
  filename = clin_txt,
  destdir = work_dir,
  getGPL = FALSE
)

clin <- Biobase::pData(geo)

colnames(clin) <- gsub("\\W", "_", colnames(clin))

data.table::fwrite(
  clin,
  file = clin_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Clinical metadata saved: ", clin_output)

# ------------------------------------------------------------
# 4. Read expression matrix
# ------------------------------------------------------------

message("Reading expression matrix...")

expr <- data.table::fread(
  text = paste(
    readLines(gzfile(expr_gz, "rt")),
    collapse = "\n"
  ),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 5. Map expression samples to GEO accessions
# ------------------------------------------------------------

sample_columns <- setdiff(colnames(expr), "Symbol_ID")

match_index <- match(sample_columns, clin$title)

if (anyNA(match_index)) {
  stop("Unmatched expression sample titles: ",
       paste(sample_columns[is.na(match_index)],
             collapse = ", "))
}

geo_ids <- clin$geo_accession[match_index]

data.table::setnames(
  expr,
  old = sample_columns,
  new = geo_ids
)

# ------------------------------------------------------------
# 6. Export expression matrix
# ------------------------------------------------------------

data.table::fwrite(
  expr,
  file = expr_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Expression matrix saved: ", expr_output)
message("Genes/features: ", nrow(expr))
message("Expression samples: ", length(geo_ids))
message("Preprocessing completed successfully.")
