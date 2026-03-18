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

## Usage

Download the files and use locally:

``` r

# Download BAM files
download.file(
  "https://zenodo.org/records/19072084/files/CD4_T_sorted_10x.bam",
  destfile = "CD4_T_sorted_10x.bam"
)

# Process with fastF
FastQDesign::fastF_bam2db(
  bam = "CD4_T_sorted_10x.bam",
  feature = "genes.tsv",
  barcode = "whitelist.txt",
  out = "output/"
)

# Load Seurat object
library(Seurat)
seurat_obj <- readRDS("CD4_T_sorted_10x.rds")
```

## Citation

If you use these datasets, please cite:

    FastQDesign Test Data (2024). Zenodo. 
    https://zenodo.org/records/19072084
    https://zenodo.org/records/19073177
