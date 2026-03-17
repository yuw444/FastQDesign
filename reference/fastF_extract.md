# Extract Custom Tags from BAM

Extract Custom Tags from BAM

## Usage

``` r
fastF_extract(bam, tag, type = 0L)
```

## Arguments

- bam:

  Path to BAM file (required)

- tag:

  Tag to extract, e.g., "CB" (required)

- type:

  Tag type: 0=string, 1=integer (default: 0)

## Value

Invisible NULL. Creates tag_summary.csv in current directory.
