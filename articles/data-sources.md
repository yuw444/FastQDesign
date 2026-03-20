# External Data Sources

Data files for FastQDesign are hosted on Zenodo for reproducibility and
easy access. Note: Some large data files (like BAM files and RDS
objects) are not included in the R package due to size limits.

## Download Links

### BAM Files (scRNA-seq alignments)

Download from: <https://zenodo.org/records/19073177>

Contains Cellranger-aligned BAM files with cell barcodes (CB tag) and
UMI (UB tag):

| File                                 | Size    |
|--------------------------------------|---------|
| `P1_AI_possorted_genome_bam.bam`     | 13.1 GB |
| `P1_AI_possorted_genome_bam.bam.bai` | 4.1 MB  |
| `P1_BM_possorted_genome_bam.bam`     | 14.4 GB |
| `P1_BM_possorted_genome_bam.bam.bai` | 4.4 MB  |

### Reference Data and Power Analysis

Download from: <https://zenodo.org/records/19072084>

Contains pre-processed Seurat objects and power analysis data:

| File | Description | Size |
|----|----|----|
| `reference_list.rds` | Reference SamplePrep output | 191.2 MB |
| `bam_downsample_list.rds` | Processed downsampled data | 77.9 MB |
| `AIBM_power.csv` | Power analysis data for experiment design | ~6 KB |

## URL Streaming (No Download Required)

BAM files can be accessed directly from Zenodo URL without downloading.

### Using Rsamtools

``` r

library(Rsamtools)

bam_url <- "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam"
bai_url <- "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam.bai"

bf <- BamFile(bam_url, index = bai_url)
header <- scanBamHeader(bf)
```

### Using fastF

fastF supports URL streaming for `extract` and `crb` subcommands:

``` r

# Extract CR tags (UMI sequences)
# Processes 182M reads in ~8 minutes from URL
FastQDesign::fastF_extract(
  bam = "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam",
  tag = "CR",
  type = 0
)
```

The `extract` subcommand streams BAM via htslib’s HTTP support. For very
large files, this is efficient as only necessary regions are downloaded.

## Local Download

``` r

# Download BAM
download.file(
  "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam",
  destfile = "P1_AI_possorted_genome_bam.bam"
)

# Download Seurat object
download.file(
  "https://zenodo.org/records/19072084/files/reference_list.rds",
  destfile = "reference_list.rds"
)

library(Seurat)
seurat_obj <- readRDS("reference_list.rds")
```

## Citation

    Wang Y (2026). FastQDesign Test Data. Zenodo.
    https://zenodo.org/records/19072084
    https://zenodo.org/records/19073177
