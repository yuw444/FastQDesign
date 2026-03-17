# Find the optimal number of cell clusters for Seurat object

Find the optimal number of cell clusters using `fpc::prediction.strenth`

## Usage

``` r
OptimalNumClusters(seu, reduction = c("pca", "umap"))
```

## Arguments

- seu:

  A seurat objet

- reduction:

  Reduction to use as input for `fpc:prediction.strength`

## Value

An object of class `predstr`, please see `fpc` for the detail
