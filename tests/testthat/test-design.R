test_that("FastQDesign works", {

  df_power_paper <- read.csv("~/Documents/Research/FastQDesign/AIBM_power.csv")

  budget <- 7500
  power_threshold <- 0.8
  flow_capacities <- c(10^7, 5 * 10^7, 2 * 10^8)
  flow_costs <- c(1000, 2000, 3000)
  library_costs <- 5000

  library(dplyr)
  library(ggplot2)
  rst_design <- FastQDesign(
    df_power = df_power_paper %>% dplyr::select(n_cells, expected_reads_per_cell, power_fastq),
    budget = budget,
    power_threshold = power_threshold,
    flow_capacities = flow_capacities,
    flow_costs = flow_costs,
    library_costs = library_costs
  )

  ggsave("./AI_BM_share_design.pdf",
         rst_design$p_share,
         width = 11,
         height = 6
  )

  p_design <- patchwork::wrap_plots(
    list(rst_design$p_ind + theme(legend.position = "none"), rst_design$p_share),
    ncol = 1
  ) + patchwork::plot_layout(guides = "collect")

  ggsave("./AI_BM_fastq_design.pdf",
         plot = p_design,
         width = 10,
         height = 10
  )

})
