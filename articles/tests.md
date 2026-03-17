# Test Suite

## Overview

FastQDesign has two test suites:

1.  **R Tests** (`tests/testthat/`) - Test R package functions
2.  **C Tests** (`tests/c/`) - Test fastF C tool

## R Tests

R tests use testthat and cover:

| Test File           | Description                 |
|---------------------|-----------------------------|
| `test-fastF.R`      | fastF subcommand wrappers   |
| `test-design.R`     | Experiment design functions |
| `test-downsample.R` | Downsampling functions      |
| `test-pseudotime.R` | Pseudotime analysis         |

### Running R Tests

``` r

devtools::test()
# Or run specific test file
testthat::test_file("tests/testthat/test-fastF.R")
```

## C Tests

C tests use Unity test framework and cover:

| Test File            | Description                     |
|----------------------|---------------------------------|
| `test_filter.c`      | FastQ filtering, whitelist tree |
| `test_utils.c`       | Sampling, sequencing utilities  |
| `test_bam2db.c`      | BAM to SQLite conversion        |
| `test_integration.c` | All fastF subcommands           |

### Running C Tests

From the build directory:

``` bash
cd build
make test_filter
./test_filter

make test_utils
./test_utils

make test_bam2db
./test_bam2db

make test_integration
./test_integration
```

Or run all tests:

``` bash
cd build
ctest --output-on-failure
```

## Test Data

Test data is located in `data/`:

- `test_R1.fastq.gz`, `test_R2.fastq.gz` - Test FastQ files
- `test.bam` - Test BAM file
- `whitelist.txt` - Cell barcode whitelist
- `genes.tsv` - Feature list
