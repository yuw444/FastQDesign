test_that("sample prep", {
  library(Seurat)

  cb03 <- readRDS("~/Documents/Research/FastQDesign/cb03.rds")

  cb03_list <- SamplePrep(
    cb03,
    condition = "orig.ident",
    cell_3d_embedding = TRUE,
    interactive = TRUE,
    min.pct = 0.2,
    logfc.threshold = 0.3,
    return.thresh = 0.05
  )

  cb03_list$Seurat$

  bam_ds_list <- readRDS("~/Documents/Research/FastQDesign/bam_downsample_list.rds")[[1]]

})
