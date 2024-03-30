test_that("interactive", {

  library(monocle3)
  library(Seurat)

  bam_ds_list <- readRDS("~/Documents/Research/FastQDesign/bam_downsample_list.rds")[[1]]

  temp <- FindPseudotime(bam_ds_list, interactive = TRUE)

  temp <- FindPseudotime(bam_ds_list)

})
