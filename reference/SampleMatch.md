# Match Downsample to Reference

Match the features of downsample to the reference, cluster membership,
differential expressed genes by cluster and condition, pseudotime, cell
lower dimension embedding

## Usage

``` r
SampleMatch(reference_list, downsample_list, ...)
```

## Arguments

- reference_list:

  The list structure same with `SamplePrep` returns

- downsample_list:

  The list structure same with `SamplePrep` returns

- ...:

  Parameters that pass to `MarkerGeneFilter`

## Value

List of matched features
