test_that("sample prep", {
  library(Seurat)

  ref_list <- readRDS("/scratch/g/chlin/Yu/FastQDesign/data/reference_list.rds")

  root_cells <- c(
    names(ref_list$dfPsedutimeCd4$pseudotime)[ref_list$dfPsedutimeCd4$root_cells],
    names(ref_list$dfPsedutimeCd8$pseudotime)[ref_list$dfPsedutimeCd8$root_cells]
  )

  ds <- DownSample(
    ref_list[[1]],
    rate_cells = 0.1,
    enable_PCR = TRUE
  )

  cb03_list <- SamplePrep(
    ref_list[[1]],
    use_default_res = TRUE,
    # condition = "orig.ident",
    cell_3d_embedding = FALSE,
    interactive = FALSE,
    min.pct = 0.2,
    logfc.threshold = 0.3,
    return.thresh = 0.05,
    verbose = TRUE
  )

  bam_ds_list <- readRDS("/scratch/g/chlin/Yu/FastQDesign/data/bam_downsample_list.rds")

  ds_list <- SamplePrep(
    bam_ds_list[[1]],
    n_clusters = 4,
    condition = "orig.ident",
    cell_3d_embedding = TRUE,
    root_cells_ref = root_cells,
    min.pct = 0.2,
    logfc.threshold = 0.3,
    return.thresh = 0.05,
    verbose = TRUE
  )

  match_list <- SampleMatch(cb03_list, ds_list)

})
