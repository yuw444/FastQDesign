# Find Pseudotime for `Seurat` Object

Find the pseudotime of a `Seurat` object based on its `UMAP`.

## Usage

``` r
FindPseudotime(
  seu,
  redo_sctransform = FALSE,
  vars_to_regress = NULL,
  use_partition = FALSE,
  close_loop = FALSE,
  learn_graph_control = NULL,
  interactive = FALSE,
  root_cells_ref = NA
)
```

## Arguments

- seu:

  A `Seurat` object

- redo_sctransform:

  Whether to redo the
  [`Seurat::SCTransform`](https://satijalab.org/seurat/reference/SCTransform.html)

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

  Whether to choose the root node interactively

- root_cells_ref:

  An vector of cell ids can be used as the root of cell trajectory

## Value

A `cell_data_set` object with pseudotime as feature
