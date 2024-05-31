#' @title FastQDesign
#' @description Design single-cell RNA experiment with raw FastQ reads
#'
#' @importFrom magrittr %>%
#' @importFrom scam scam
#' @importFrom purrr map map2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr filter arrange mutate desc slice_head
#' @importFrom ggplot2 ggplot geom_point scale_fill_gradientn scale_color_manual
#'      scale_y_continuous scale_x_continuous labs scale_shape_manual theme_bw
#'      geom_hline geom_vline geom_text guides theme
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom scales label_number cut_short_scale
#' @param df_power A data frame with power information,
#'    \code{N}(filtered cell number), the cell numbers after filtration in the data analysis(\code{Seurat}) pipeline
#'    \code{R}(reads required to get per filtered cell), the total number of FastQ reads divide by the filtered cell numbers
#'            the total number may use \code{fastF} to get for each sample
#'    \code{power}(0-1) defined by the weighted average of the Adjusted Rand Index on matched cluster membership,
#'          the Jaccard Index on cluster DE genes and condition DE genes,
#'          the correlation index on the matched pseudotime(optional),
#'          the correlation index on the matched 3d cell embedding(optional)
#'
#' @param budget A numeric of budget
#' @param power_threshold A numeric of power threshold
#' @param reads_valid_rate The percentage of FastQ reads that is valid when converted to UMIs in the FastQ reference
#'  \code{./filtered_feature_bc_matrix/barcodes.tsv.gz} after the alignment
#' @param flow_capacities A vector of flow capacities
#' @param flow_costs A vector of flow costs
#' @param library_costs A vector of library costs
#' @param cell_increment A numeric of cell increment for the design, default is \code{floor(max(df_power$N) / 50) * 5}
#' @param read_increment A numeric of read increment for the design, default is \code{floor(max(df_power$R) / 100) * 10}
#' @return A list with power information under constraints and ggplots
#'
#' @export
#'
# Generated from function body. Editing this file has no effect.
FastQDesign <-
  function (df_power,
            budget,
            power_threshold,
            reads_valid_rate,
            flowcell_capacities,
            flow_costs,
            library_costs,
            cell_increment = NA,
            read_increment = NA)
  {
    if (budget < min(library_costs + flow_costs)) {
      stop(
        strwrap(
          "\033[31mError: Budget is smaller than\n      minimum cost(`flow_costs` + `library_costs`), please\n      increase it!\033[0m"
        )
      )
    }
    df_power <- as.data.frame(df_power)
    colnames(df_power) <- c("N", "R", "power")
    if (power_threshold < min(df_power$power)) {
      stop(
        strwrap(
          "\033[31mError: `power_threshold` is bigger\n      than the maximum power in `df_power`, please\n      decrease it!\033[0m"
        )
      )
    }
    original_capacities <- flowcell_capacities
    flowcell_capacities <- flowcell_capacities * reads_valid_rate
    if (is.na(cell_increment)) {
      cell_increment <- floor(max(df_power$N) / 50) * 5
    }
    if (is.na(read_increment)) {
      read_increment <- floor(max(df_power$R) / 100) * 10
    }
    N <- seq(cell_increment,
             floor(max(df_power$N) / cell_increment) *
               cell_increment,
             cell_increment)
    R_shared_base <-
      seq(read_increment,
          floor(max(df_power$R) / read_increment) *
            read_increment,
          read_increment)
    max_N <- max(N)
    max_R <- max(R_shared_base)
    if (max(flowcell_capacities) < min(df_power$N * df_power$R)) {
      warning(
        strwrap(
          "\033[31m The maximum of `flowcell_capacities` is smaller\n      than the minimum of `N*R` in `df_power`, the result\n      may be biased!\033[0m"
        )
      )
    }
    model_scam <- scam::scam(power ~ s(N, k = 10, bs = "mpi") +
                               s(R, k = 10, bs = "mpi"), data = df_power)
    df_flow_list <- list()
    for (i in seq_along(flowcell_capacities)) {
      R <- flowcell_capacities[i] / N
      df_flow_list[[i]] <-
        data.frame(
          N = N,
          R = R,
          power = predict(model_scam,
                          newdata = data.frame(N = N, R = R))
        ) %>% dplyr::filter(N <
                              max_N, R < max_R) %>% dplyr::mutate(flow_cell = i,
                                                                  cost = library_costs + flow_costs[i])
    }
    df_flow <-
      do.call(rbind, df_flow_list) %>% dplyr::mutate(flow_cell = factor(flow_cell)) %>%
      dplyr::mutate(
        feasible = ifelse(cost <= budget, TRUE, FALSE),
        feasible = factor(
          feasible,
          levels = c(TRUE, FALSE),
          labels = c(TRUE, FALSE)
        )
      )
    temp <- function(y) {
      function(x) {
        y / x
      }
    }
    capacties_afforable <-
      (budget - library_costs) / min(flow_costs / original_capacities) *
      reads_valid_rate
    curves_flow_shared <-
      purrr::map(c(flowcell_capacities, capacties_afforable),
                 ~ temp(.x))
    converted_capacities <-
      (scales::label_number(scale_cut = scales::cut_short_scale()))(original_capacities)
    levels(df_flow$flow_cell) <- converted_capacities
    ind_available_design <- df_flow %>% dplyr::filter(power >=
                                                        power_threshold, cost <= budget)
    if (nrow(ind_available_design) == 0) {
      warning(
        strwrap(
          "\033[33m Missing optimal indepedent design with current budget\n        and power restriction, please increase budget or lower\n        power restriction!\033[0m"
        )
      )
    }
    ind_optimal_power <-
      ind_available_design %>% dplyr::arrange(dplyr::desc(power)) %>%
      dplyr::slice_head() %>% dplyr::mutate(type = "similarity")
    ind_optimal_cost <-
      ind_available_design %>% dplyr::arrange(cost,
                                              dplyr::desc(power)) %>%
      dplyr::slice_head() %>% dplyr::mutate(type = "cost")
    ind_optimal <- rbind(ind_optimal_power, ind_optimal_cost)
    colors_ind <-
      RColorBrewer::brewer.pal(n = length(flowcell_capacities) +
                                 3, name = "Set1")[seq_len(1 + length(flowcell_capacities))]
    names(colors_ind) <- c(converted_capacities, "Budget")
    type_values <- c(similarity = 8, cost = 9)
    p_design_ind <- ggplot2::ggplot(df_flow, ggplot2::aes(N,
                                                          R)) +
      ggplot2::geom_point(shape = 21,
                          size = 3,
                          ggplot2::aes(fill = feasible)) +
      ggplot2::scale_fill_manual(
        values = c("grey30", "grey70"),
        breaks = c(TRUE, FALSE),
        drop = FALSE
      ) +
      purrr::map2(
        curves_flow_shared,
        c(converted_capacities, "Budget"),
        ~ ggplot2::stat_function(fun = .x,
                                 ggplot2::aes(color = .y))
      ) + ggplot2::scale_color_manual(
        values = colors_ind,
        breaks = c(converted_capacities, "Budget"),
        name = "Flow Cell Capacity"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(suffix = " K",
                                      scale = 0.001),
        limits = c(0, max(df_flow$R) + 10)
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::label_number(suffix = " K",
                                      scale = 0.001),
        limits = c(0, max(df_flow$N) + 10)
      ) +
      ggplot2::labs(x = "Target number of cells", y = "Target reads per cell", fill = "Feasibility") +
      ggplot2::geom_point(
        data = ind_optimal[!is.na(ind_optimal$N),],
        ggplot2::aes(N, R, shape = type),
        size = 4,
        color = "purple"
      ) +
      ggplot2::scale_shape_manual(values = type_values, name = "Optimal Design") +
      ggplot2::theme_bw()
    p_power_cost_ind <- ggplot2::ggplot(df_flow, ggplot2::aes(cost,
                                                              power)) + ggplot2::geom_point(shape = 21,
                                                                                            size = 3,
                                                                                            ggplot2::aes(fill = feasible)) +
      ggplot2::scale_fill_manual(
        values = c("grey30", "grey70"),
        breaks = c(TRUE, FALSE),
        drop = FALSE
      ) +
      ggplot2::labs(x = "Cost",
                    y = "Similarity", fill = "Feasibility") +
      ggplot2::geom_hline(yintercept = power_threshold,
                          col = "red") +
      ggplot2::geom_vline(xintercept = budget,
                          col = "blue") + ggplot2::geom_text(
                            x = budget,
                            y = min(df_flow$power) +
                              diff(range(df_flow$power)) /
                              10,
                            label = "Budget",
                            col = "blue"
                          ) +
      ggplot2::geom_text(
        x = max(df_flow$cost, na.rm = TRUE) -
          diff(range(df_flow$cost) / 5),
        y = power_threshold,
        label = "Threshold",
        col = "red"
      ) + ggplot2::geom_point(
        data = ind_optimal[!is.na(ind_optimal$N),],
        ggplot2::aes(cost, power, shape = type),
        size = 4,
        color = "purple"
      ) + ggplot2::scale_shape_manual(values = type_values,
                                      name = "Optimal Design") +
      ggplot2::guides(color = "none") +
      ggplot2::theme_bw()
    p_ind <- patchwork::wrap_plots(
      p_design_ind + p_power_cost_ind +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    )
    grid_list <- list(N = N, R = R_shared_base)
    df_power_shared <- expand.grid(grid_list)
    df_power_shared$power <-
      predict(model_scam, newdata = df_power_shared)
    shared_flow_cells_index <-
      which.min(flow_costs / flowcell_capacities)
    flow_cell_used <- converted_capacities[shared_flow_cells_index]
    df_power_shared$cost <-
      library_costs + min(flow_costs / flowcell_capacities) *
      df_power_shared$R * df_power_shared$N
    df_power_shared$flow_cell <- flow_cell_used
    df_power_shared$feasible <- ifelse(df_power_shared$cost <=
                                         budget, TRUE, FALSE)
    df_power_shared$feasible <- factor(
      df_power_shared$feasible,
      levels = c(TRUE, FALSE),
      labels = c(TRUE, FALSE)
    )
    share_available_design <-
      df_power_shared %>% dplyr::filter(power >=
                                          power_threshold, cost <= budget)
    if (nrow(share_available_design) == 0) {
      warning(
        strwrap(
          "\033[33m Missing optimal share design with current budget\n        and power restriction, please increase budget or lower\n        power restriction!\033[0m"
        )
      )
    }
    share_optimal_power <-
      share_available_design %>% dplyr::arrange(dplyr::desc(power)) %>%
      dplyr::slice_head() %>% dplyr::mutate(type = "similarity")
    share_optimal_cost <-
      share_available_design %>%
      dplyr::arrange(cost,
                     dplyr::desc(power)) %>%
      dplyr::slice_head() %>%
      dplyr::mutate(type = "cost")
    share_optimal <- rbind(share_optimal_power, share_optimal_cost)
    p_design_share <-
      ggplot2::ggplot(df_power_shared, ggplot2::aes(N,
                                                    R)) +
      ggplot2::geom_point(shape = 21,
                          size = 3,
                          ggplot2::aes(fill = feasible)) +
      ggplot2::scale_fill_manual(
        values = c("grey30", "grey70"),
        breaks = c(TRUE, FALSE),
        drop = FALSE
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(suffix = " K",
                                      scale = 0.001),
        limits = c(0, max(df_power_shared$R) +
                     10)
      ) + ggplot2::scale_x_continuous(
        labels = scales::label_number(suffix = " K",
                                      scale = 0.001),
        limits = c(0, max(df_power_shared$N) +
                     10)
      ) + purrr::map2(
        curves_flow_shared,
        c(converted_capacities,
          "Budget"),
        ~ ggplot2::stat_function(fun = .x, ggplot2::aes(color = .y))
      ) +
      ggplot2::scale_color_manual(
        values = colors_ind,
        breaks = c(converted_capacities,
                   "Budget"),
        name = "Flow Cell Capacity"
      ) + ggplot2::geom_point(
        data = share_optimal,
        ggplot2::aes(N, R, shape = type),
        size = 4,
        color = "purple"
      ) +
      ggplot2::scale_shape_manual(values = type_values, name = "Optimal Design") +
      ggplot2::labs(
        x = "Target number of cells",
        y = "Target reads per cell",
        shape = "Optimal Design",
        fill = "Feasibility"
      ) + ggplot2::theme_bw()
    p_power_cost_share <-
      ggplot2::ggplot(df_power_shared, ggplot2::aes(cost,
                                                    power)) +
      ggplot2::geom_point(shape = 21,
                          size = 3,
                          ggplot2::aes(fill = feasible)) +
      ggplot2::scale_fill_manual(
        values = c("grey30", "grey70"),
        breaks = c(TRUE, FALSE),
        drop = FALSE
      ) +
      ggplot2::labs(x = "Cost",
                    y = "Similarity", fill = "Feasibility") +
      ggplot2::geom_hline(yintercept = power_threshold,
                          col = "red") +
      ggplot2::geom_vline(xintercept = budget,
                          col = "blue") +
      ggplot2::geom_text(
        x = budget,
        y = min(df_flow$power) +
          diff(range(df_flow$power)) /
          10,
        label = "Budget",
        col = "blue"
      ) +
      ggplot2::geom_text(
        x = max(df_flow$cost, na.rm = TRUE) -
          diff(range(df_flow$cost) / 5),
        y = power_threshold,
        label = "Threshold",
        col = "red"
      ) + ggplot2::geom_point(
        data = share_optimal,
        ggplot2::aes(cost, power, shape = type),
        size = 4,
        color = "purple"
      ) + ggplot2::scale_shape_manual(values = type_values,
                                      name = "Optimal Design") +
      ggplot2::guides(color = "none") +
      ggplot2::theme_bw()
    p_share <-
      patchwork::wrap_plots(
        p_design_share + p_power_cost_share +
          patchwork::plot_layout(guides = "collect") &
          ggplot2::theme(legend.position = "bottom")
      )
    output <-
      list(
        ind_design = df_flow,
        share_design = df_power_shared,
        p_ind = p_ind,
        p_share = p_share,
        p_design_ind = p_design_ind,
        p_design_share = p_design_share,
        p_power_cost_ind = p_power_cost_ind,
        p_power_cost_share = p_power_cost_share,
        ind_available_design = ind_available_design,
        ind_optimal_cost = ind_optimal_cost,
        ind_optimal_power = ind_optimal_power,
        share_available_design = share_available_design,
        share_optimal_cost = share_optimal_cost,
        share_optimal_power = share_optimal_power
      )
    return(output)
  }
