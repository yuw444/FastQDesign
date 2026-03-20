# Generating Power Metadata with fastF

This vignette demonstrates how to generate `df_power` metadata for
experiment design using the `FastQDesign-fastF` tool and FastQDesign
framework. The `df_power` data frame contains power/similarity metrics
across a grid of cell proportions and read depths.

## Overview

The workflow consists of:

1.  **Prepare reference data** - BAM files from AI and BM samples
2.  **Generate downsampled UMI matrices** - Using
    [`fastF_bam2db()`](https://yuw444.github.io/FastQDesign/reference/fastF_bam2db.md)
    across a parameter grid
3.  **Compare with reference** - Using
    [`SamplePrep()`](https://yuw444.github.io/FastQDesign/reference/SamplePrep.md)
    and
    [`SampleMatch()`](https://yuw444.github.io/FastQDesign/reference/SampleMatch.md)
4.  **Aggregate results** - Calculate mean similarity across replicates

## Required Inputs

- BAM files with aligned reads (`possorted_genome_bam.bam`)
- Feature annotation (`features.tsv.gz`)
- Cell barcode list (`barcodes.tsv.gz`)
- Reference UMI matrix for comparison

## Data Source

The AI and BM BAM files used in this vignette are available from Zenodo:

**BAM Files:** <https://zenodo.org/records/19073177>

| Sample | File                             | Size    |
|--------|----------------------------------|---------|
| AI     | `P1_AI_possorted_genome_bam.bam` | 13.1 GB |
| BM     | `P1_BM_possorted_genome_bam.bam` | 14.4 GB |

**Pre-processed Seurat Objects:** <https://zenodo.org/records/19072084>

| File                      | Description                 | Size     |
|---------------------------|-----------------------------|----------|
| `reference_list.rds`      | Reference SamplePrep output | 191.2 MB |
| `bam_downsample_list.rds` | Processed downsampled data  | 77.9 MB  |

### Downloading BAM Files

``` bash
# Download AI BAM
wget https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam
wget https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam.bai

# Download BM BAM
wget https://zenodo.org/records/19073177/files/P1_BM_possorted_genome_bam.bam
wget https://zenodo.org/records/19073177/files/P1_BM_possorted_genome_bam.bam.bai
```

Or use R to download:

``` r

# Download AI BAM
download.file(
  "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam",
  destfile = "P1_AI_possorted_genome_bam.bam"
)
download.file(
  "https://zenodo.org/records/19073177/files/P1_AI_possorted_genome_bam.bam.bai",
  destfile = "P1_AI_possorted_genome_bam.bam.bai"
)

# Download BM BAM
download.file(
  "https://zenodo.org/records/19073177/files/P1_BM_possorted_genome_bam.bam",
  destfile = "P1_BM_possorted_genome_bam.bam"
)
download.file(
  "https://zenodo.org/records/19073177/files/P1_BM_possorted_genome_bam.bai",
  destfile = "P1_BM_possorted_genome_bam.bai"
)
```

### Generating Feature and Barcode Files

The BAM files contain Cell Ranger-aligned reads with CB (cell barcode)
and UB (UMI) tags. After downloading the Seurat objects from Zenodo, you
can extract the feature and barcode files:

``` r

library(Seurat)

# Load reference Seurat object
ref_list <- readRDS("reference_list.rds")

# Extract barcodes from reference
writeLines(colnames(ref_list$Seurat), "ai_barcodes.tsv")
writeLines(colnames(ref_list$Seurat), "bm_barcodes.tsv")

# For features, use Cell Ranger reference (e.g., mm10 or GRCh38)
# Or extract from the BAM header using samtools:
# samtools view -H P1_AI_possorted_genome_bam.bam | grep "^@SQ"
```

For more details on data preparation, see the [External Data
Sources](https://yuw444.github.io/FastQDesign/data-sources/index.md)
vignette.

## Step 1: Define Configuration

``` r

library(FastQDesign)
library(Seurat)
library(dplyr)
library(future)
plan("multisession", workers = 4)

# Paths to reference data (update these paths after downloading from Zenodo)
# BAM files: https://zenodo.org/records/19073177
AI_BAM <- "./P1_AI_possorted_genome_bam.bam"
BM_BAM <- "./P1_BM_possorted_genome_bam.bam"

# Feature and barcode files - extract from Cell Ranger output directory
# or download from the BAM file's companion files
AI_FEATURES <- "./P1_AI_filtered_feature_bc_matrix/features.tsv.gz"
BM_FEATURES <- "./P1_BM_filtered_feature_bc_matrix/features.tsv.gz"
AI_BARCODES <- "./P1_AI_filtered_feature_bc_matrix/barcodes.tsv.gz"
BM_BARCODES <- "./P1_BM_filtered_feature_bc_matrix/barcodes.tsv.gz"

# Alternative: Load pre-processed reference from Zenodo
# https://zenodo.org/records/19072084
# ref_list <- readRDS("reference_list.rds")

# Output directory
OUT_DIR <- "./power_analysis"

# Grid parameters (reduced to 5 replicates per grid point)
cell_rates <- seq(0.1, 1.0, by = 0.1)  # 10 cell proportions
read_rates <- seq(0.1, 1.0, by = 0.1)  # 10 read depths
seeds <- 926:930  # 5 replicates (was 10 in original: 926:935)

# Reference data (full UMI matrices)
# If using Cell Ranger output, point to filtered_feature_bc_matrix directories
REF_AI <- "./P1_AI_filtered_feature_bc_matrix"
REF_BM <- "./P1_BM_filtered_feature_bc_matrix"
```

## Step 2: Prepare Reference

Load and prepare the reference UMI matrices:

``` r

# Load reference UMI matrices
AI_data <- Read10X(REF_AI)
BM_data <- Read10X(REF_BM)

# Create Seurat objects
AI <- CreateSeuratObject(counts = AI_data, project = "ai", min.cells = 1, min.features = 1)
BM <- CreateSeuratObject(counts = BM_data, project = "bm", min.cells = 1, min.features = 1)

# Merge samples
AIBM <- merge(AI, BM, add.cell.ids = c("ai", "bm"), project = "AIBM")
AIBM[["percent.mt"]] <- PercentageFeatureSet(AIBM, pattern = "^mt-")

# Prepare reference for comparison
ref_list <- SamplePrep(
  AIBM,
  n_clusters = 4,
  pca_used = 1:30,
  use_default_res = FALSE,
  condition = "orig.ident",
  vars_to_regress = "percent.mt",
  cell_3d_embedding = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.3,
  return.thresh = 0.05,
  verbose = TRUE
)
```

## Step 3: Generate Downsampled Matrices with fastF

Loop through the grid and generate downsampled UMI matrices:

``` r

# Create output directories
dir.create(file.path(OUT_DIR, "fastF"), recursive = TRUE, showWarnings = FALSE)

# Grid dimensions
n_cells <- length(cell_rates)  # 10
n_reads <- length(read_rates)   # 10
n_reps <- length(seeds)        # 5

total_jobs <- n_cells * n_reads * n_reps
cat("Total grid points:", n_cells * n_reads, "\n")
cat("Replicates per grid:", n_reps, "\n")
cat("Total jobs:", total_jobs, "\n")

# Function to generate one downsampled matrix
generate_downsample <- function(seed, cell_rate, read_rate, sample_name, 
                                 bam, feature, barcode, out_dir) {
  
  out_path <- file.path(out_dir, paste0("seed", seed, "/N", cell_rate, "/R", read_rate, "/", sample_name))
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  # Run fastF bam2db
  FastQDesign::fastF_bam2db(
    bam = bam,
    feature = feature,
    barcode = barcode,
    out = out_path,
    dbname = ":memory:",
    cell = cell_rate,
    depth = read_rate,
    seed = seed
  )
  
  return(out_path)
}

# Iterate through grid
for (seed in seeds) {
  for (cell_rate in cell_rates) {
    for (read_rate in read_rates) {
      
      cat("Processing: seed=", seed, "cell=", cell_rate, "read=", read_rate, "\n")
      
      # Generate AI downsampled matrix
      generate_downsample(seed, cell_rate, read_rate, "AI", 
                          AI_BAM, AI_FEATURES, AI_BARCODES, OUT_DIR)
      
      # Generate BM downsampled matrix  
      generate_downsample(seed, cell_rate, read_rate, "BM",
                          BM_BAM, BM_FEATURES, BM_BARCODES, OUT_DIR)
    }
  }
}
```

## Step 4: Compare Downsampled to Reference

Load each downsampled matrix, compare with reference, and compute
similarity metrics:

``` r

# Function to compute similarity metrics
compute_similarity <- function(input_path, ref_list) {
  
  # Read downsampled matrices
  AI_ds <- Read10X(file.path(input_path, "AI"))
  BM_ds <- Read10X(file.path(input_path, "BM"))
  
  # Create Seurat objects
  AI <- CreateSeuratObject(counts = AI_ds, project = "ai", min.cells = 1, min.features = 1)
  BM <- CreateSeuratObject(counts = BM_ds, project = "bm", min.cells = 1, min.features = 1)
  
  # Merge
  AIBM_ds <- merge(AI, BM, add.cell.ids = c("ai", "bm"), project = "AIBM")
  AIBM_ds[["percent.mt"]] <- PercentageFeatureSet(AIBM_ds, pattern = "^mt-")
  
  # Subset to reference cells (if downsampled below reference)
  AIBM_ds <- subset(AIBM_ds, cells = Cells(ref_list$Seurat))
  
  # Prepare downsampled
  ds_list <- SamplePrep(
    AIBM_ds,
    n_clusters = 4,
    pca_used = 1:30,
    use_default_res = FALSE,
    condition = "orig.ident",
    vars_to_regress = "percent.mt",
    cell_3d_embedding = TRUE,
    min.pct = 0.2,
    logfc.threshold = 0.3,
    return.thresh = 0.05,
    verbose = FALSE
  )
  
  # Match with reference
  match_list <- SampleMatch(ref_list, ds_list, avg_log2FC_abs = 0.3, p_val_adj = 0.05)
  
  # Compute metrics
  ARI <- mclust::adjustedRandIndex(
    match_list$cluster_match$seurat_clusters.x,
    match_list$cluster_match$seurat_clusters.y
  )
  
  JI_cluster <- JaccardIndex(
    match_list$genes_cluster$ref,
    match_list$genes_cluster$ds
  )
  
  JI_condition <- JaccardIndex(
    match_list$genes_condition$ref,
    match_list$genes_condition$ds
  )
  
  Kendall_cell <- pcaPP::cor.fk(
    match_list$cluster_match$pseudotime.x,
    match_list$cluster_match$pseudotime.y
  )
  
  # Kendall for 3D embedding
  Kendall_embed <- pcaPP::cor.fk(
    dist(match_list$cluster_match[, c("umap_1.x", "umap_2.x", "umap_3.x")]),
    dist(match_list$cluster_match[, c("umap_1.y", "umap_2.y", "umap_3.y")])
  )
  
  # Mean similarity (power)
  power <- mean(c(ARI, JI_cluster, JI_condition, Kendall_cell, Kendall_embed))
  
  return(power)
}

# Collect all results
results <- list()

for (seed in seeds) {
  for (cell_rate in cell_rates) {
    for (read_rate in read_rates) {
      input_path <- file.path(OUT_DIR, paste0("seed", seed, "/N", cell_rate, "/R", read_rate))
      
      tryCatch({
        power <- compute_similarity(input_path, ref_list)
        
        results[[length(results) + 1]] <- data.frame(
          seed = seed,
          cell_rate = cell_rate,
          read_rate = read_rate,
          power = power,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
      })
    }
  }
}

df_results <- bind_rows(results)
```

## Step 5: Aggregate Replicates and Format df_power

Average across replicates and convert to the format expected by
[`FastQDesign()`](https://yuw444.github.io/FastQDesign/reference/FastQDesign.md):

``` r

# Aggregate: mean power across replicates per grid point
df_power_agg <- df_results %>%
  group_by(cell_rate, read_rate) %>%
  summarize(
    mean_power = mean(power, na.rm = TRUE),
    sd_power = sd(power, na.rm = TRUE),
    .groups = "drop"
  )

# Get reference cell and read counts
ref_n_cells <- ncol(ref_list$Seurat)
ref_total_reads <- sum(ref_list$Seurat$nCount_RNA)
ref_reads_per_cell <- ref_total_reads / ref_n_cells

# Convert rates to absolute values
df_power <- df_power_agg %>%
  mutate(
    n_cells = round(cell_rate * ref_n_cells),
    expected_reads_per_cell = round(read_rate * ref_reads_per_cell)
  ) %>%
  select(n_cells, expected_reads_per_cell, mean_power) %>%
  rename(power_fastq = mean_power)

# Filter out grid points with insufficient cells
df_power <- df_power %>% filter(n_cells > 0, expected_reads_per_cell > 0)

# Save results
write.csv(df_power, file.path(OUT_DIR, "df_power.csv"), row.names = FALSE)

cat("\nGenerated df_power with", nrow(df_power), "grid points\n")
head(df_power)
```

## Step 6: Use df_power for Experiment Design

Once you have the `df_power` data frame, you can use it with
[`FastQDesign()`](https://yuw444.github.io/FastQDesign/reference/FastQDesign.md):

``` r

# Define experiment parameters
budget <- 7500
power_threshold <- 0.7
flowcell_capacities <- c(10^7, 5 * 10^7, 2 * 10^8)
flowcell_costs <- c(1000, 2000, 3000)
library_costs <- 5000

# Run experiment design
rst <- FastQDesign(
  df_power = df_power,
  budget = budget,
  power_threshold = power_threshold,
  reads_valid_rate = 0.9,
  flowcell_capacities = flowcell_capacities,
  flowcell_costs = flowcell_costs,
  library_costs = library_costs
)

# View results
print(rst$ind_optimal)
print(rst$share_optimal)
```

## Output Format

The `df_power` data frame should have the following columns:

| Column                    | Description                     |
|---------------------------|---------------------------------|
| `n_cells`                 | Number of cells after filtering |
| `expected_reads_per_cell` | Total reads / filtered cells    |
| `power_fastq`             | Mean similarity score (0-1)     |

## Batch Processing with SLURM

For large grids, submit as a SLURM array job:

``` bash
#!/bin/bash
#SBATCH --job-name=power_grid
#SBATCH --array=1-500
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=12:00:00

# Calculate indices
SEEDS=(926 927 928 929 930)
CELLS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
READS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

idx=$((SLURM_ARRAY_TASK_ID - 1))
seed_idx=$((idx % 5))
cell_idx=$(( (idx / 5) % 10 ))
read_idx=$(( idx / 50 ))

SEED=${SEEDS[$seed_idx]}
CELL=${CELLS[$cell_idx]}
READ=${READS[$read_idx]}

# Run fastF for both samples
fastF bam2db -b AI.bam -f features.tsv -a barcodes.tsv \
    -c $CELL -r $READ -s $SEED -o output/seed${SEED}/N${CELL}/R${READ}/AI

fastF bam2db -b BM.bam -f features.tsv -a barcodes.tsv \
    -c $CELL -r $READ -s $SEED -o output/seed${SEED}/N${CELL}/R${READ}/BM
```

## Notes

- The grid uses rates (0-1) which are converted to absolute values based
  on reference
- 5 replicates (seeds) reduce computation while maintaining statistical
  power
- Similarity metrics include: ARI, Jaccard Index (cluster genes),
  Jaccard Index (condition genes), Kendall correlation (pseudotime),
  Kendall correlation (3D embedding)
- Power/similarity of 1.0 indicates perfect match to reference
