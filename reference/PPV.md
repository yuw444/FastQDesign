# Positive Predictive Value(PPV)

Calculate the PPV index between two sets

## Usage

``` r
PPV(a, b)
```

## Arguments

- a:

  an vector or list, predicted value

- b:

  an vector or list, true value

## Value

scalar, PPV index

## Examples

``` r
set.seed(926)
a <- sample(1:100, 50)
b <- sample(1:100, 50)
PPV(a,b)
#> [1] 0.52
```
