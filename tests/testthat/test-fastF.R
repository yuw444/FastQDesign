# Test fastF subcommands from R using FastQDesign package wrappers
# This file tests all 5 fastF CLI subcommands via the R package:
# - fastF_freq: Find cell barcode whitelist and frequencies from FastQ
# - fastF_filter: Filter FastQ files using whitelist and read depth
# - fastF_crb: Extract CR/CB tags from BAM and summarize to TSV
# - fastF_bam2db: Convert BAM to SQLite UMI matrix
# - fastF_extract: Extract custom tags from BAM to CSV

DATA_PATH <- "/scratch/g/chlin/Yu/FastQDesign/data"
BUILD_PATH <- "/scratch/g/chlin/Yu/FastQDesign/build/testthat_output"

test_that("fastF_freq - find cell barcode whitelist", {
  # Skip if binary not found (will be auto-located by package)
  
  out_dir <- file.path(BUILD_PATH, "freq")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Use package wrapper function
  FastQDesign::fastF_freq(
    R1 = file.path(DATA_PATH, "test_R1.fastq.gz"),
    out = out_dir,
    len = 16L
  )
  
  # Check output file exists
  whitelist_file <- file.path(out_dir, "whitelist.txt")
  expect_true(file.exists(whitelist_file), 
              paste("Expected whitelist.txt at", whitelist_file))
  
  # Verify whitelist has content
  whitelist_lines <- readLines(whitelist_file)
  expect_gt(length(whitelist_lines), 0)
})

test_that("fastF_filter - filter FastQ files using whitelist", {
  out_dir <- file.path(BUILD_PATH, "filter")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Use package wrapper function
  FastQDesign::fastF_filter(
    R1 = file.path(DATA_PATH, "test_R1.fastq.gz"),
    R2 = file.path(DATA_PATH, "test_R2.fastq.gz"),
    whitelist = file.path(DATA_PATH, "whitelist.txt"),
    out = out_dir,
    len = 16L,
    seed = 42L,
    rate = 0.1
  )
  
  # Check output files exist
  expect_true(file.exists(file.path(out_dir, "R1.fastq.gz")),
              "R1.fastq.gz should be created")
  expect_true(file.exists(file.path(out_dir, "R2.fastq.gz")),
              "R2.fastq.gz should be created")
})

test_that("fastF_crb - extract CR/CB tags from BAM", {
  out_file <- file.path(BUILD_PATH, "crb_output.tsv.gz")
  
  # Use package wrapper function
  FastQDesign::fastF_crb(
    bam = file.path(DATA_PATH, "test.bam"),
    out = out_file
  )
  
  # Check output file exists
  expect_true(file.exists(out_file), "crb output should be created")
})

test_that("fastF_bam2db - convert BAM to SQLite UMI matrix", {
  out_dir <- file.path(BUILD_PATH, "bam2db")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Use package wrapper function
  FastQDesign::fastF_bam2db(
    bam = file.path(DATA_PATH, "test.bam"),
    feature = file.path(DATA_PATH, "features.tsv.gz"),
    barcode = file.path(DATA_PATH, "barcodes.tsv.gz"),
    out = out_dir
  )
  
  # Check output files exist (matrix market format)
  expect_true(file.exists(file.path(out_dir, "matrix.mtx.gz")),
              "matrix.mtx.gz should be created")
  expect_true(file.exists(file.path(out_dir, "barcodes.tsv.gz")),
              "barcodes.tsv.gz should be created")
  expect_true(file.exists(file.path(out_dir, "features.tsv.gz")),
              "features.tsv.gz should be created")
})

test_that("fastF_extract - extract custom tags from BAM", {
  out_dir <- file.path(BUILD_PATH, "extract")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Save current directory and change to output dir
  # (extract outputs to current working directory)
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(out_dir)
  
  # Use package wrapper function
  FastQDesign::fastF_extract(
    bam = file.path(DATA_PATH, "test.bam"),
    tag = "CB"
  )
  
  # Output file is tag_summary.csv in the current directory
  out_file <- file.path(out_dir, "tag_summary.csv")
  expect_true(file.exists(out_file), "extract output should be created")
  
  # Verify it's not empty
  lines <- readLines(out_file)
  expect_gt(length(lines), 1)
})
