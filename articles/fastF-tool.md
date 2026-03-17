# Using fastF Tool

The `fastF` tool is a C-based command-line utility for scRNA-seq
experiment design. It provides several subcommands for processing FastQ
and BAM files.

## Installation

The fastF binary is automatically installed with the FastQDesign
package:

``` r

devtools::install_github("yuw444/FastQDesign")
```

## fastF Subcommands

### 1. Find Cell Barcode Whitelist

Extract cell barcode whitelist from FastQ files:

``` r

# Find whitelist from R1 FastQ file
FastQDesign::fastF_freq(
  R1 = "test_R1.fastq.gz",
  out = "output_dir",
  len = 16L
)
```

This creates `whitelist.txt` in the output directory.

### 2. Filter FastQ Files

Filter FastQ files using a whitelist and downsampling:

``` r

FastQDesign::fastF_filter(
  R1 = "test_R1.fastq.gz",
  R2 = "test_R2.fastq.gz",
  whitelist = "whitelist.txt",
  out = "filtered_dir",
  len = 16L,
  seed = 42L,
  rate = 0.1
)
```

### 3. Extract CR/CB Tags from BAM

Summarize cell barcode and UMI information:

``` r

FastQDesign::fastF_crb(
  bam = "test.bam",
  out = "output.tsv.gz"
)
```

### 4. Convert BAM to Matrix

Convert BAM to sparse matrix format:

``` r

FastQDesign::fastF_bam2db(
  bam = "test.bam",
  feature = "genes.tsv",
  barcode = "whitelist.txt",
  out = "matrix_dir"
)
```

### 5. Extract Custom Tags

Extract custom tags from BAM files:

``` r

FastQDesign::fastF_extract(
  bam = "test.bam",
  tag = "CB"
)
```

## Running from Command Line

You can also run fastF directly from the command line:

``` bash
# Find whitelist
fastF freq -R input_R1.fastq.gz -l 16 -o output_dir/

# Filter FastQ
fastF filter -R input_R1.fastq.gz -r input_R2.fastq.gz -w whitelist.txt -o filtered/

# Extract CR/CB tags
fastF crb -b input.bam -o output.tsv.gz
```
