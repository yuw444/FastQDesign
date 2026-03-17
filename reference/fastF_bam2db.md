# Convert BAM to SQLite UMI Matrix

Convert BAM to SQLite UMI Matrix

## Usage

``` r
fastF_bam2db(
  bam,
  feature,
  barcode,
  out = ".",
  dbname = ":memory:",
  cell = 1,
  depth = 1,
  seed = 926L,
  umicopies = FALSE
)
```

## Arguments

- bam:

  Path to BAM file (required)

- feature:

  Path to feature list file (required)

- barcode:

  Path to barcode list file (required)

- out:

  Output directory (default: ".")

- dbname:

  Database name (default: ":memory:")

- cell:

  Rate of cell barcode (default: 1.0)

- depth:

  Rate of depth (default: 1.0)

- seed:

  Random seed (default: 926)

- umicopies:

  Whether to store umi.tsv.gz (default: FALSE)

## Value

Invisible NULL. Creates SQLite DB in output directory.
