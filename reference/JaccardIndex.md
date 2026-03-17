# Jaccard Index

Calculate the Jaccard index between two sets

## Usage

``` r
JaccardIndex(a, b)
```

## Arguments

- a:

  an vector or list

- b:

  an vector or list

## Value

scalar, jaccard index

## Examples

``` r
set.seed(926)
a <- sample(1:100, 50)
b <- sample(1:100, 50)
JaccardIndex(a,b)
#> [1] 0.351
```
