# Sample Preparation and Analysis

This vignette demonstrates sample preparation, downsampling, and
pseudotime analysis functions.

## Data Source

Reference data for sample preparation is hosted on Zenodo:

**Download from:** <https://zenodo.org/records/19072084>

``` r

# Download reference data
download.file(
  "https://zenodo.org/records/19072084/files/reference_list.rds",
  destfile = "reference_list.rds"
)
```

## SamplePrep

Prepare reference data for experiment design:

``` r

# Load reference data (after downloading from Zenodo)
ref_list <- readRDS("reference_list.rds")

# Prepare sample (Seurat object)
cb03_list <- SamplePrep(
  ref_list[[1]],
  use_default_res = TRUE,
  cell_3d_embedding = FALSE,
  interactive = FALSE,
  min.pct = 0.2,
  logfc.threshold = 0.3,
  return.thresh = 0.05,
  verbose = TRUE
)
```

### Parameters

| Parameter           | Description                                  |
|---------------------|----------------------------------------------|
| `object`            | Seurat object or BAM file list               |
| `use_default_res`   | Use default clustering resolution            |
| `n_clusters`        | Number of clusters (if not using default)    |
| `condition`         | Condition column for differential expression |
| `cell_3d_embedding` | Compute 3D embedding                         |
| `min.pct`           | Minimum percentage for markers               |
| `logfc.threshold`   | Log fold change threshold                    |

## DownSample

Downsample cells to simulate different sequencing depths:

``` r

# Downsample reference
ds <- DownSample(
  ref_list[[1]],
  rate_cells = 0.1,
  enable_PCR = TRUE
)
```

### Parameters

| Parameter    | Description                     |
|--------------|---------------------------------|
| `object`     | Seurat object                   |
| `rate_cells` | Fraction of cells to keep       |
| `enable_PCR` | Enable PCR duplicate simulation |

## FindPseudotime

Perform pseudotime analysis using Monocle3:

``` r

# Perform pseudotime analysis
pt_result <- FindPseudotime(
  cb03_list,
  root_cells = NULL,  # or specify root cell names
  verbose = TRUE
)
```

The function returns a list with: - **$`cds**: Monocle3 CellDataSet
- **`$pseudotime**: Pseudotime values - **\$plot**: Pseudotime
visualization

## SampleMatch

Compare experimental samples to reference:

``` r

# Compare samples
match_result <- SampleMatch(
  query = query_seurat,
  reference = reference_list,
  method = "seurat"
)
```
