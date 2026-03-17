# Root Node Selection for `cell_data_set`

Select the root node for `cell_data_set` based on the given
`root_cells_ref`

## Usage

``` r
RootNodeSelect(cds, root_cells_ref = NA)
```

## Arguments

- cds:

  A `cell_data_set` object from `monocle3`

- root_cells_ref:

  An vector of cell ids can be used as the root of cell trajectory

## Value

A character, rootnode id
