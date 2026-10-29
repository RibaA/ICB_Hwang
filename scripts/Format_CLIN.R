#!/usr/bin/env Rscript

# ============================================================
# Clinical metadata harmonization: GSE136961 (Hwang)
# ============================================================
#
# Input:
#   data/download/CLIN.txt
#   data/annot/curation_tissue.csv
#
# Output:
#   data/curation/CLIN.csv (semicolon-delimited)
#
# Purpose:
#   Add standardized ICB clinical variables while preserving
#   all original GEO metadata columns.
#
# Notes:
#   - Histology is encoded using OncoTree histology codes.
#   - Survival times are converted from days to months.
#   - The source recurrence indicator is stored as 'pfs'.
#   - The cohort received nivolumab OR pembrolizumab (anti-PD-1);
#     individual treatment assignments are not available here.
#   - 'drug_name' lists both study-reported treatment options;
#     it does not mean that each patient received both drugs.
#   - Get_Response() may infer response from PFS rather than
#     reproduce the study's RECIST-based DCB/NDB classification.
#   - Missing clinical values are left as NA.
# ============================================================

suppressPackageStartupMessages({
  library(stringr)
  library(tibble)
})

# 1. Directories ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_dir  <- if (length(args) >= 1) args[1] else "data/download"
output_dir <- if (length(args) >= 2) args[2] else "data/curation"
annot_dir  <- if (length(args) >= 3) args[3] else "data/annot"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Shared harmonization functions -----------------------------
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

# 3. Load original GEO metadata ---------------------------------
clin_raw <- read.csv(
  file.path(input_dir, "CLIN.txt"),
  stringsAsFactors = FALSE,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# 4. Add standardized columns without dropping GEO fields ------
selected_cols <- c(
  "patient", "sex", "primary", "recist", "drug_type",
  "t.os", "t.pfs", "histo", "stage", "age", "os",
  "pfs", "dna", "dna_info", "rna", "rna_info",
  "response.other.info", "response"
)

clin <- clin_raw
new_cols <- setdiff(selected_cols, colnames(clin))
clin[new_cols] <- NA

# 5. Patient, tissue, treatment and expression ------------------
clin$patient <- clin_raw$geo_accession
clin$sex <- ifelse(
  tolower(trimws(clin_raw$Sex_ch1)) == "female", "F",
  ifelse(tolower(trimws(clin_raw$Sex_ch1)) == "male", "M", NA_character_)
)
clin$primary <- "Lung"
clin$drug_type <- "PD-1"
clin$rna <-  "rnaseq"
clin$rna_info <-  "tpm"

# Harmonize histology using descriptive OncoTree names
histo_map <- c(
  "adenocarcinoma" = "Lung Adenocarcinoma",
  "adenocarcinoma-squamous cell carcinoma" = "Lung Adenosquamous Carcinoma",
  "Large cell neuroendocarine carcinoma" = "Lung Large Cell Neuroendocrine Carcinoma",
  "squamous cell carcinoma" = "Lung Squamous Cell Carcinoma"
)

clin$histo <- unname(histo_map[clin_raw$tumor_type_ch1])

# 6. Survival information ---------------------------------------
# GEO text example:
# recurrence:1, progression-free-survival_days:144, live,
# overall-survival_days:234
survival <- stringr::str_split(clin_raw$survival_ch1, ",")

clin$pfs <- as.numeric(vapply(survival, function(s) {
  if (length(s) < 1) return(NA_character_)
  stringr::str_extract(s[1], "\\d+")
}, character(1)))

clin$t.pfs <- as.numeric(vapply(survival, function(s) {
  if (length(s) < 2) return(NA_character_)
  stringr::str_extract(s[2], "\\d+")
}, character(1))) / 30.5

clin$os <- vapply(survival, function(s) {
  if (length(s) < 3) return(NA_real_)
  status <- tolower(stringr::str_trim(s[3]))
  if (status == "live") return(0)
  if (status == "dead") return(1)
  NA_real_
}, numeric(1))

clin$t.os <- as.numeric(vapply(survival, function(s) {
  if (length(s) < 4) return(NA_character_)
  stringr::str_extract(s[4], "\\d+")
}, character(1))) / 30.5

# 7. Shared response and tissue annotations ---------------------
clin$response <- Get_Response(data = clin)

annotation_tissue <- read.csv(
  file.path(annot_dir, "curation_tissue.csv"),
  stringsAsFactors = FALSE
)
clin <- annotate_tissue(
  clin = clin,
  study = "Hwang",
  annotation_tissue = annotation_tissue,
  check_histo = FALSE
)

# Annotating Tissue Data
annotation_tissue <- read.csv("https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv")
clin <- annotate_tissue(clin=clin, study='Hwang', annotation_tissue= annotation_tissue, check_histo=FALSE)

# Set treatmentid based on curation_drug.csv file.
#annotation_drug <- read.csv("https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_drug.csv")
#clin <- add_column(clin, treatmentid=clin$drug_type, .after='tissueid')
clin <- tibble::add_column(
  clin,
  treatmentid = "Nivolumab|Pembrolizumab",
  .after = "tissueid"
)


# 8. Place standardized variables first, keep all GEO fields ----
clin <- clin[, c(selected_cols,
                 setdiff(colnames(clin), selected_cols)),
             drop = FALSE]

# 9. Export -----------------------------------------------------
output_file <- file.path(output_dir, "CLIN.csv")
write.table(
  clin,
  file = output_file,
  quote = FALSE,
  sep = ";",
  col.names = TRUE,
  row.names = FALSE,
  na = "NA"
)
message("Harmonized clinical metadata saved: ", output_file)
message("Clinical samples: ", nrow(clin))
message("Clinical harmonization completed successfully.")
