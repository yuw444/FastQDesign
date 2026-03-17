# Down Sample a Seurat object

Down sample a seurat object with given \`rate_cells\` and \`rate_umis\`

## Usage

``` r
DownSample(
  seu,
  rate_cells = 0.3,
  rate_umis = 0.3,
  seed = 926,
  enable_PCR = FALSE,
  nb_size = 2,
  nb_prob = 0.2
)
```

## Arguments

- seu:

  A Seurat object

- rate_cells:

  The proportion of cell to sample

- rate_umis:

  The proportion of UMIs to sample

- seed:

  A seed for random number generation

- enable_PCR:

  Whether to consider the duplication of UMIs, default is \`TRUE\`, use
  negative binomial distribution to simulate the number of copies of
  each UMI

- nb_size:

  If \`enable_PCR=TRUE\`, it is used in
  [rnbinom](https://rdrr.io/r/stats/NegBinomial.html)

- nb_prob:

  If \`enable_PCR=TRUE\`, it is used in
  [rnbinom](https://rdrr.io/r/stats/NegBinomial.html)
