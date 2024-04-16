test_that("sample prep", {
  library(Seurat)

  ref_list <- readRDS("~/Documents/Research/FastQDesign/reference_list.rds")

  ds <- DownSample(
    ref_list[[1]],
    rate_cells = 0.1,
    enable_PCR = TRUE
  )

  cb03_list <- SamplePrep(
    ref_list[[1]],
    condition = "orig.ident",
    cell_3d_embedding = TRUE,
    interactive = TRUE,
    min.pct = 0.2,
    logfc.threshold = 0.3,
    return.thresh = 0.05,
    verbose = TRUE
  )

  bam_ds_list <- readRDS("~/Documents/Research/FastQDesign/bam_downsample_list.rds")

  ds_list <- SamplePrep(
    bam_ds_list[[1]],
    n_clusters = 4,
    condition = "orig.ident",
    cell_3d_embedding = TRUE,
    root_cells_ref = Cells(cb03_list[[1]])[cb03_list[[1]]$root_cells],
    min.pct = 0.2,
    logfc.threshold = 0.3,
    return.thresh = 0.05,
    verbose = TRUE
  )

  match_list <- SampleMatch(cb03_list, ds_list)

})
