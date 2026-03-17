# Find Cell Barcode Whitelist from FastQ

Find Cell Barcode Whitelist from FastQ

## Usage

``` r
fastF_freq(R1, out = ".", len = 16L, umi = 10L)
```

## Arguments

- R1:

  Path to R1 fastq files (required)

- out:

  Output directory for whitelist (default: ".")

- len:

  Length of cell barcode (default: 16)

- umi:

  Length of UMI (default: 10)

## Value

Invisible NULL. Creates whitelist.txt in output directory.
