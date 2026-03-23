# FastQDesign

<!-- badges: start -->
<!-- badges: end -->

The goal of FastQDesign is to guide the investigator in designing a single-cell RNA sequencing(scRNA-seq) experiment. We aim to shift the focus back to raw FastQ reads other than the Unique Molecular Identifier (UMI) matrix when considering the scRNA-seq experiment. The [paper](https://www.nature.com/articles/s42003-025-07938-8) is now open-access, PMID: [40175506](https://pubmed.ncbi.nlm.nih.gov/40175506/). 

## Installation

### Install R Package

You can install the development version of FastQDesign from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("yuw444/FastQDesign")
```

### Install fastF Tool (Optional)

The `fastF` C tool is automatically installed with the R package. For manual installation or building from source:

```bash
# Clone the repository
git clone https://github.com/yuw444/FastQDesign.git
cd FastQDesign

# Build fastF (requires htslib, zlib, libbz2, openssl, sqlite3)
mkdir build && cd build
cmake ..
make

# Add to PATH
export PATH=$PATH:/path/to/FastQDesign/build/fastF
```

### System Dependencies

- **R** (>= 3.5.0)
- **htslib** (>= 1.10)
- **zlib**, **libbz2**, **openssl**, **sqlite3**

On Ubuntu/Debian:
```bash
sudo apt-get install libhts-dev libbz2-dev libssl-dev libsodium-dev libsqlite3-dev
```

On macOS:
```bash
brew install htslib bzip2 openssl sqlite
```

## fastF

[fastF](https://github.com/yuw444/fastF) is a submodel of the FastQDesign framework, it is written in C for efficiency and broad compatibility. It helps generate pseudo-design datasets from the FastQ reference. 

## Prepare the Reference and Subsample

This is a basic example which shows you how to prepare the subsamples and reference for the comparison:

```
library(FastQDesign)
library(Seurat)

# Load reference data from Zenodo (https://zenodo.org/records/19072084)
download.file(
  "https://zenodo.org/records/19072084/files/reference_list.rds",
  destfile = "reference_list.rds"
)
aibm <- readRDS("reference_list.rds")[[1]]

ref_list <- SamplePrep(
  aibm,
  condition = "orig.ident",
  cell_3d_embedding = TRUE,
  interactive = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.3,
  return.thresh = 0.05,
  verbose = TRUE
)

# Download downsampled data from Zenodo
download.file(
  "https://zenodo.org/records/19072084/files/bam_downsample_list.rds",
  destfile = "bam_downsample_list.rds"
)
fastq_ds <- readRDS("bam_downsample_list.rds")[[1]]

ds_list <- SamplePrep(
  fastq_ds,
  n_clusters = 4,
  condition = "orig.ident",
  cell_3d_embedding = TRUE,
  root_cells_ref = Cells(aibm)[aibm$root_cells],
  min.pct = 0.2,
  logfc.threshold = 0.3,
  return.thresh = 0.05,
  verbose = TRUE
)

match_list <- SampleMatch(cb03_list, ds_list)

```

In the above example, `aibm` is the reference Seurat object, `fastq_ds` is generated with the FastQ downsample by [fastF](https://github.com/yuw444/fastF).


## Generate One Dot in the Similarity Surface

```
## the similiarity from one subsample

ARI <- mclust::adjustedRandIndex(match_list$seurat_cluster.x, match_list$seurat_cluster.y)

JaccardCluster <- FastQDesign:::JaccardIndex(match_list$genes_cluster$ref, match_list$genes_cluster$ds)
JaccardCondition <- FastQDesign:::JaccardIndex(match_list$genes_condition$ref, match_list$genes_condition$ds)

Kendall <- cor(match_list$cluster_match$pseudotime.x, match_list$cluster_match$pseudotime.y, method = "kendall")

similarity <- mean(c(ARI, JaccardCluster, JaccardCondition, Kendall))
```

## Design

```
# Load power analysis data from package
df_power_paper <- read.csv(system.file("data", "AIBM_power.csv", package = "FastQDesign"))

budget <- 7500
power_threshold <- 0.8
flow_capacities <- c(10^7, 5 * 10^7, 2 * 10^8)
flow_costs <- c(1000, 2000, 3000)
library_costs <- 5000

library(dplyr)
library(ggplot2)
rst_design <- FastQDesign(
  df_power = df_power_paper %>% dplyr::select(n_cells, expected_reads_per_cell, power_fastq),
  budget = budget,
  power_threshold = power_threshold,
  flow_capacities = flow_capacities,
  flow_costs = flow_costs,
  library_costs = library_costs
)
```
