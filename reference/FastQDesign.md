# FastQDesign

Design single-cell RNA experiment with raw FastQ reads

## Usage

``` r
FastQDesign(
  df_power,
  budget,
  power_threshold,
  reads_valid_rate,
  flowcell_capacities,
  flowcell_costs,
  library_costs,
  cell_increment = NA,
  read_increment = NA
)
```

## Arguments

- df_power:

  A data frame with power information, `N`(filtered cell number), the
  cell numbers after filtration in the data analysis(`Seurat`) pipeline
  `R`(reads required to get per filtered cell), the total number of
  FastQ reads divide by the filtered cell numbers the total number may
  use `fastF` to get for each sample `power`(0-1) defined by the
  weighted average of the Adjusted Rand Index on matched cluster
  membership, the Jaccard Index on cluster DE genes and condition DE
  genes, the correlation index on the matched pseudotime(optional), the
  correlation index on the matched 3d cell embedding(optional)

- budget:

  A numeric of budget

- power_threshold:

  A numeric of power threshold

- reads_valid_rate:

  The percentage of FastQ reads that is valid when converted to UMIs in
  the FastQ reference `./filtered_feature_bc_matrix/barcodes.tsv.gz`
  after the alignment

- flowcell_capacities:

  A vector of flow capacities

- flowcell_costs:

  A vector of flow costs

- library_costs:

  A vector of library costs

- cell_increment:

  A numeric of cell increment for the design, default is
  `floor(max(df_power$N) / 50) * 5`

- read_increment:

  A numeric of read increment for the design, default is
  `floor(max(df_power$R) / 100) * 10`

## Value

A list with power information under constraints and ggplots
