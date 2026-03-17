# Extract CR/CB Tags from BAM

Extract CR/CB Tags from BAM

## Usage

``` r
fastF_crb(bam, out = "crb_output.tsv")
```

## Arguments

- bam:

  Path to BAM file (required)

- out:

  Output TSV file (default: "crb_output.tsv")

## Value

Invisible NULL. Creates TSV file with CB-\>CR-\>count mapping.
