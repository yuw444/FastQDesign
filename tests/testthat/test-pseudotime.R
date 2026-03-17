test_that("interactive", {
  library(monocle3)
  library(Seurat)

  ref_list <- readRDS(
    "/scratch/g/chlin/Yu/FastQDesign/data/reference_list.rds"
  )

  bam_ds <- readRDS(
    "/scratch/g/chlin/Yu/FastQDesign/data/bam_downsample_list.rds"
  )[[1]]

  root_cells <- c(
    names(ref_list$dfPsedutimeCd4$pseudotime)[
      ref_list$dfPsedutimeCd4$root_cells
    ],
    names(ref_list$dfPsedutimeCd8$pseudotime)[
      ref_list$dfPsedutimeCd8$root_cells
    ]
  )

  temp <- FindPseudotime(
    bam_ds,
    root_cells_ref = root_cells
  )

  # temp <- FindPseudotime(bam_ds)

  cb03_cds <- FindPseudotime(
    ref_list[[1]],
    interactive = FALSE,
    root_cells_ref = root_cells,
    use_partition = TRUE
  )

  plot_cells(cb03_cds, color_cells_by = "pseudotime")
})
