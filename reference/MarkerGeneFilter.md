# DE Markers Filters

Filter the maker genes with series of filters

## Usage

``` r
MarkerGeneFilter(
  df_MarkerGene,
  pct_1 = 0.2,
  pct_2 = 0.2,
  p_val_adj_ = 0.05,
  avg_log2FC_abs = 0
)
```

## Arguments

- df_MarkerGene:

  a \`data.fram\` generated from \`FindMarkers\`

- pct_1:

  A minimum cutoff for \`pct.1\`

- pct_2:

  A minimum cutoff for \`pct.2\`

- p_val_adj\_:

  A maximum cutoff for \`p_val_adj\`

- avg_log2FC_abs:

  A minimum cutoff for \`abs(avg_log2FC)\`

## Value

An vector of genes that pass through the filters
