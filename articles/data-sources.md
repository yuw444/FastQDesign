# External Data Sources

Test data for FastQDesign is hosted on Zenodo for reproducibility and
easy access.

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

### Seurat Objects (processed data)

Download from: <https://zenodo.org/records/19072084>

Contains pre-processed Seurat objects:

| File                      | Size     |
|---------------------------|----------|
| `bam_downsample_list.rds` | 77.9 MB  |
| `reference_list.rds`      | 191.2 MB |

## URL Streaming (No Download Required)

BAM files can be accessed directly from Zenodo URL. Note: Zenodo BAM
files are large (13-14 GB each), so streaming may be slow.

### Using Rsamtools

``` r

library(Rsamtools)

bam_url <- "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam"
bai_url <- "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam.bai"

bf <- BamFile(bam_url, index = bai_url)
header <- scanBamHeader(bf)
```

### Using fastF

fastF supports URL streaming but downloading is recommended for large
BAM files:

``` r

# Works but slow for 13GB files
FastQDesign::fastF_crb(
  bam = "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam",
  out = "output.tsv.gz"
)
```

**Recommendation:** For large BAM files, download first rather than
stream.

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
