# External Data Sources

Test data for FastQDesign is hosted on Zenodo for reproducibility and
easy access.

## Download Links

### BAM Files (scRNA-seq alignments)

Download from: <https://zenodo.org/records/19072084>

Contains BAM files with cell barcodes (CB tag) and UMI (UB tag) for
testing: - `CD4_T_sorted_10x.bam` - `CD8_T_sorted_10x.bam` -
`PBMC_10x.bam`

### Seurat Objects (processed data)

Download from: <https://zenodo.org/records/19073177>

Contains pre-processed Seurat objects with: - Normalized counts -
PCA/UMAP embeddings - Cluster annotations - Cell metadata

## Direct URL Streaming

BAM files can be accessed directly from URL without downloading,
provided: 1. The remote server supports HTTP range requests 2. Both
`.bam` and `.bai` (index) files are accessible via URL

### Using Rsamtools

``` r

library(Rsamtools)

# Remote BAM with remote index
bam_url <- "https://zenodo.org/records/19072084/files/CD4_T_sorted_10x.bam"
bai_url <- "https://zenodo.org/records/19072084/files/CD4_T_sorted_10x.bam.bai"

# Create BamFile with remote access
bf <- BamFile(bam_url, index = bai_url)

# Scan reads
scanBam(bf)
```

### Using fastF from URL

``` r

# Download BAM files
FastQDesign::fastF_bam2db(
  bam = "https://zenodo.org/records/19072084/files/CD4_T_sorted_10x.bam",
  feature = "genes.tsv",
  barcode = "whitelist.txt",
  out = "output/"
)
```

## Local Download

Alternatively, download files first:

``` r

# Download BAM files
download.file(
  "https://zenodo.org/records/19072084/files/CD4_T_sorted_10x.bam",
  destfile = "CD4_T_sorted_10x.bam"
)
download.file(
  "https://zenodo.org/records/19072084/files/CD4_T_sorted_10x.bam.bai",
  destfile = "CD4_T_sorted_10x.bam.bai"
)

# Download Seurat object
download.file(
  "https://zenodo.org/records/19073177/files/CD4_T_sorted_10x.rds",
  destfile = "CD4_T_sorted_10x.rds"
)

# Load in R
library(Seurat)
seurat_obj <- readRDS("CD4_T_sorted_10x.rds")
```

## Citation

If you use these datasets, please cite:

    FastQDesign Test Data (2024). Zenodo. 
    https://zenodo.org/records/19072084
    https://zenodo.org/records/19073177
