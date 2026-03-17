# Conditional Shannon Entropy

Calculate Conditional Shannon Entropy, weighted by its proportion on the
condition

## Usage

``` r
ConditionalShannonEntropy(condition, feature)
```

## Arguments

- condition:

  an vector, the condition

- feature:

  an vector, the interested feature

## Value

a scalar

## Examples

``` r
set.seed(926)
condition <- sample(1:100, 50)
feature <- sample(1:100, 50)
ConditionalShannonEntropy(condition,feature)
#> [1] 0
```
