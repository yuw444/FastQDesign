# Filter FastQ Files Using Whitelist

Filter FastQ Files Using Whitelist

## Usage

``` r
fastF_filter(
  R1,
  R2 = NULL,
  I1 = NULL,
  out = ".",
  whitelist = NULL,
  len = 16L,
  seed = 926L,
  rate = 0,
  allcells = FALSE
)
```

## Arguments

- R1:

  Path to R1 fastq files (required)

- R2:

  Path to R2 fastq files (optional)

- I1:

  Path to I1 fastq files (optional)

- out:

  Output directory (default: ".")

- whitelist:

  Path to whitelist file (required unless allcells=TRUE)

- len:

  Length of cell barcode (default: 16)

- seed:

  Random seed (default: 926)

- rate:

  Rate of reads to keep (default: 0)

- allcells:

  Keep all reads with cell barcode (default: FALSE)

## Value

Invisible NULL. Creates I1.fastq.gz, R1.fastq.gz, R2.fastq.gz in output.
