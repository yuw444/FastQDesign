# Find Markers By Condition

Find Differential Expressed(DE) markers by Condition

## Usage

``` r
FindAllMarkersByCondition(seu, condition, ...)
```

## Arguments

- seu:

  A Seurat object

- condition:

  A character name in the \`meta.data\` of \`seu\`

- ...:

  The argument could pass to
  [`FindMarkers`](https://satijalab.org/seurat/reference/FindMarkers.html)

## Value

A data frame includes DE markers
