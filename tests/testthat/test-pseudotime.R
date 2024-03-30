test_that("interactive", {

  library(monocle3)
  library(Seurat)

  bam_ds <- readRDS("~/Documents/Research/FastQDesign/bam_downsample_list.rds")[[1]]

  temp <- FindPseudotime(bam_ds, interactive = TRUE)

  temp <- FindPseudotime(bam_ds)

  cb03_list <- readRDS("~/Documents/Research/FastQDesign/reference_list.rds")

  cb03_cds <- FindPseudotime(cb03_list[[1]], interactive = TRUE, use_partition = TRUE)

  plot_cells(cb03_cds, color_cells_by = "pseudotime")

})
