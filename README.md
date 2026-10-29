# ICB_Hwang: Immunotherapy Cohort Data Curation

Reproducible data curation and harmonization of the **Hwang immunotherapy cohort** ([GEO GSE136961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136961)) using **R** and **Snakemake**. The workflow generates harmonized clinical metadata, expression data, assay-availability information, and a `MultiAssayExperiment` (MAE) object using **GENCODE v19 (hg19)**.

## Repository structure

```text
ICB_Hwang/
├── README.md
├── Snakefile
├── scripts/
│   ├── format_downloaded_data.R
│   ├── Format_CLIN.R
│   ├── Format_cased_sequenced.R
│   ├── Format_EXPR.R
│   └── create_MAE.R
└── data/
    ├── download/             # GEO inputs and preprocessed files
    ├── annot/                # Tissue and GENCODE v19 annotations
    └── curation/             # Harmonized outputs and MAE
```

## Data curation workflow

The `Snakefile` runs the following steps based on file dependencies:

| Step | Snakemake rule | Script / task | Output |
| --- | --- | --- | --- |
| 1 | `download_data` | Download GEO clinical and TPM files | `data/download/*.gz` |
| 2 | `format_downloaded_data` | Extract GEO metadata and map expression sample IDs to GEO accessions | `CLIN.txt`, `EXP_TPM.tsv` |
| 3 | `download_annotation` | Obtain the shared tissue annotation | `curation_tissue.csv` |
| 4 | `format_clin` | Harmonize clinical fields while preserving original GEO metadata | `CLIN.csv` |
| 5 | `format_cased_sequenced` | Record SNV, CNA, and expression availability | `cased_sequenced.csv` |
| 6 | `format_expr` | Harmonize gene symbols, resolve duplicates, and log-transform TPM | `EXPR.csv` |
| 7 | `get_MultiAssayExp` | Build an MAE using GENCODE v19 gene annotations | `ICB_Hwang.rds` |

### Scripts

- **`format_downloaded_data.R`** — prepares clinical metadata and the TPM matrix from GEO downloads.
- **`Format_CLIN.R`** — adds common ICB variables (including histology, treatment, survival, and response) without dropping original GEO columns; unavailable values remain `NA`.
- **`Format_cased_sequenced.R`** — records molecular assay availability per sample.
- **`Format_EXPR.R`** — extracts gene symbols, retains one row per gene, and applies `log2(TPM + 0.001)`.
- **`create_MAE.R`** — combines clinical data and annotated expression measurements into a `MultiAssayExperiment`.

Further study-specific harmonization details can be documented in comments within each script or in `docs/`.

## Run the workflow

Requires **Snakemake 9.x**, **R** and the R/Bioconductor packages used in `scripts/`, plus `wget` for downloads. Run in Linux or WSL from the repository root.

Place `Gencode.v19.annotation.RData` (containing `features_gene`) in `data/annot/` before building the MAE; this file is supplied locally.

```bash
snakemake -n             # Preview jobs
snakemake --cores 1      # Run the workflow
```

To generate only the MAE from existing curated inputs:

```bash
snakemake --cores 1 data/curation/ICB_Hwang.rds
```

## Outputs

```text
data/curation/
├── CLIN.csv
├── cased_sequenced.csv
├── EXPR.csv
└── ICB_Hwang.rds
```

The completed Hwang MAE contains an expression `RangedSummarizedExperiment` with **381 annotated genes and 21 samples** (current run). Original GEO metadata are retained in the clinical data.

## References

- [Hwang cohort — GEO GSE136961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136961)
- [ICB_Common — shared harmonization functions](https://github.com/BHKLAB-Pachyderm/ICB_Common)
- [Snakemake documentation](https://snakemake.readthedocs.io/)
