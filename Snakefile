
# ============================================================
# ICB_Hwang (GSE136961): Data curation and harmonization
# ============================================================
#
# Workflow:
#   1. Download GEO data
#   2. Preprocess clinical and expression data
#   3. Download tissue annotation
#   4. Harmonize clinical metadata
#   5. Generate molecular data availability
#   6. Harmonize gene expression
#   7. Generate MultiAssayExperiment (MAE)
#
# Genome annotation: Gencode v19 (hg19)
# Snakemake: Version 9.x
#
# Run:
#   snakemake -n
#   snakemake --cores 1
# ============================================================

DOWNLOAD = "data/download"
CURATION = "data/curation"
ANNOT = "data/annot"


# ------------------------------------------------------------
# 1. Final outputs
# ------------------------------------------------------------

rule all:
    input:
        f"{CURATION}/CLIN.csv",
        f"{CURATION}/EXPR.csv",
        f"{CURATION}/cased_sequenced.csv",
        f"{CURATION}/ICB_Hwang.rds"


# ------------------------------------------------------------
# 2. Download GEO data
# ------------------------------------------------------------

rule download_data:
    output:
        expr=f"{DOWNLOAD}/GSE136961_TPM.tsv.gz",
        clin=f"{DOWNLOAD}/GSE136961_series_matrix.txt.gz"
    shell:
        """
        mkdir -p {DOWNLOAD}

        wget -O {output.expr:q} \
          'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/suppl/GSE136961_TPM.tsv.gz'

        wget -O {output.clin:q} \
          'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/matrix/GSE136961_series_matrix.txt.gz'
        """


# ------------------------------------------------------------
# 3. Preprocess downloaded data
# ------------------------------------------------------------

rule format_downloaded_data:
    input:
        expr=f"{DOWNLOAD}/GSE136961_TPM.tsv.gz",
        clin=f"{DOWNLOAD}/GSE136961_series_matrix.txt.gz",
        script="scripts/format_downloaded_data.R"
    output:
        clin=f"{DOWNLOAD}/CLIN.txt",
        expr=f"{DOWNLOAD}/EXP_TPM.tsv"
    shell:
        """
        Rscript {input.script:q} {DOWNLOAD}
        """


# ------------------------------------------------------------
# 4. Download tissue annotation
# ------------------------------------------------------------

rule download_annotation:
    output:
        tissue=f"{ANNOT}/curation_tissue.csv"
    shell:
        """
        mkdir -p {ANNOT}

        wget -O {output.tissue:q} \
          'https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/data/curation_tissue.csv'
        """


# ------------------------------------------------------------
# 5. Harmonize clinical metadata
# ------------------------------------------------------------

rule format_clin:
    input:
        clin=f"{DOWNLOAD}/CLIN.txt",
        tissue=f"{ANNOT}/curation_tissue.csv",
        script="scripts/format_CLIN.R"
    output:
        f"{CURATION}/CLIN.csv"
    shell:
        """
        mkdir -p {CURATION}

        Rscript {input.script:q} \
          {DOWNLOAD} \
          {CURATION} \
          {ANNOT}
        """


# ------------------------------------------------------------
# 6. Generate molecular data availability
# ------------------------------------------------------------

rule format_cased_sequenced:
    input:
        clin=f"{DOWNLOAD}/CLIN.txt",
        script="scripts/format_cased_sequenced.R"
    output:
        f"{CURATION}/cased_sequenced.csv"
    shell:
        """
        mkdir -p {CURATION}

        Rscript {input.script:q} \
          {DOWNLOAD} \
          {CURATION}
        """


# ------------------------------------------------------------
# 7. Harmonize gene expression
# ------------------------------------------------------------

rule format_expr:
    input:
        expr=f"{DOWNLOAD}/EXP_TPM.tsv",
        case=f"{CURATION}/cased_sequenced.csv",
        script="scripts/format_EXPR.R"
    output:
        f"{CURATION}/EXPR.csv"
    shell:
        """
        mkdir -p {CURATION}

        Rscript {input.script:q} \
          {DOWNLOAD} \
          {CURATION}
        """


# ------------------------------------------------------------
# 8. Generate MultiAssayExperiment (MAE)
# ------------------------------------------------------------

rule get_MultiAssayExp:
    input:
        clin=f"{CURATION}/CLIN.csv",
        expr=f"{CURATION}/EXPR.csv",
        case=f"{CURATION}/cased_sequenced.csv",
        annotation=f"{ANNOT}/Gencode.v19.annotation.RData",
        script="scripts/create_MAE.R"
    output:
        f"{CURATION}/ICB_Hwang.rds"
    shell:
        """
        Rscript {input.script:q} \
          {CURATION} \
          {input.annotation:q} \
          {output:q}
        """
