
#!/usr/bin/env Rscript

# ============================================================
# Sequencing data availability: GSE136961 (Hwang)
# ============================================================
#
# Input:
#   data/download/CLIN.txt
#
# Output:
#   data/curation/cased_sequenced.csv (semicolon-delimited)
#
# Purpose:
#   Generate a patient-level matrix indicating the availability
#   of molecular data for each sample in the Hwang cohort.
#
# Variables:
#   patient : GEO sample accession (GSM identifier)
#   snv     : Somatic mutation data availability (0/1)
#   cna     : Copy number alteration data availability (0/1)
#   expr    : Gene expression data availability (0/1)
#
# Notes:
#   - All samples have RNA expression data (expr = 1).
#   - SNV and CNA data are not provided (snv = 0, cna = 0).
#   - 1 indicates available; 0 indicates unavailable.
#   - Output follows the common ICB curation format.
# ============================================================

# ------------------------------------------------------------
# 1. Define directories
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

input_dir  <- if (length(args) >= 1) args[1] else "data/download"
output_dir <- if (length(args) >= 2) args[2] else "data/curation"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read clinical metadata
# ------------------------------------------------------------

clin <- read.csv(
  file.path(input_dir, "CLIN.txt"),
  stringsAsFactors = FALSE,
  sep = "\t",
  header = TRUE
)

# ------------------------------------------------------------
# 3. Extract patient identifiers
# ------------------------------------------------------------

patient <- sort(unique(clin$geo_accession))

# ------------------------------------------------------------
# 4. Define molecular data availability
# ------------------------------------------------------------

case <- data.frame(
  patient = patient,
  snv = 0,
  cna = 0,
  expr = 1,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 5. Export sequencing availability
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  "cased_sequenced.csv"
)

write.table(
  case,
  file = output_file,
  quote = FALSE,
  sep = ";",
  col.names = TRUE,
  row.names = FALSE
)

message("Sequencing availability saved: ", output_file)
message("Total samples: ", nrow(case))
message("Completed successfully.")
