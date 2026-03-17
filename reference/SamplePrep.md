# Prepare a sample

Prepare a sample for the comparison

## Usage

``` r
SamplePrep(
  seu,
  cluster_id = NA,
  n_clusters = NA,
  use_default_res = TRUE,
  condition = NA,
  cell_3d_embedding = FALSE,
  pca_used = 1:30,
  vars_to_regress = NULL,
  use_partition = FALSE,
  close_loop = FALSE,
  learn_graph_control = NULL,
  interactive = FALSE,
  root_cells_ref = NA,
  verbose = FALSE,
  ...
)
```

## Arguments

- seu:

  A `Seurat` object

- cluster_id:

  A A character, the condition colname in `Seurat@meta.data`

- n_clusters:

  A scalar, the desired number of clusters

- use_default_res:

  Whether to use the default resolution in
  [`Seurat::FindClusters`](https://satijalab.org/seurat/reference/FindClusters.html)

- condition:

  A character, the condition colname in `Seurat@meta.data`

- cell_3d_embedding:

  Calculate 3d embedding for the input `seu`

- pca_used:

  Used in
  [`Seurat::FindNeighbors`](https://satijalab.org/seurat/reference/FindNeighbors.html)
  and
  [`Seurat::RunUMAP`](https://satijalab.org/seurat/reference/RunUMAP.html)

- vars_to_regress:

  A character string, used in
  [`Seurat::SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)

- use_partition:

  Whether to `use_partition` when
  [`monocle3::learn_graph`](https://rdrr.io/pkg/monocle3/man/learn_graph.html)

- close_loop:

  Logic, a parameter in
  [`monocle3::learn_graph`](https://rdrr.io/pkg/monocle3/man/learn_graph.html)

- learn_graph_control:

  List, a parameter in
  [`monocle3::learn_graph`](https://rdrr.io/pkg/monocle3/man/learn_graph.html)

- interactive:

  Logic, whether to select root cell interactively when considering cell
  fate

- root_cells_ref:

  An vector of cell ids can be used as the root of cell trajectory

## Value

list Seurat with extra meta.data features, data.frames of DE genes by
cluster and condition
