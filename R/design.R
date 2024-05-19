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
#' @param cells_used_rate The percentage of Cells passed quality control in the
#'  \code{./filtered_feature_bc_matrix/barcodes.tsv.gz} after the alignment
#' @param flow_capacities A vector of flow capacities
#' @param flow_costs A vector of flow costs
#' @param library_costs A vector of library costs
#' @return A list with power information under constraints and ggplots
#'
#' @export
#'
FastQDesign <- function(df_power,
                        budget,
                        power_threshold,
                        reads_valid_rate,
                        cells_used_rate,
                        flow_capacities,
                        flow_costs,
                        library_costs) {

  # check budget constraints
  if (budget < min(library_costs + flow_costs)) {
    stop(
      strwrap(
        "\033[31mError: Budget is smaller than
      minimum cost(`flow_costs` + `library_costs`), please
      increase it!\033[0m"
      )
    )
  }

  # check power constraints
  df_power <- as.data.frame(df_power)
  colnames(df_power) <- c("N", "R", "power")
  if (power_threshold < min(df_power$power)) {
    stop(
      strwrap(
        "\033[31mError: `power_threshold` is bigger
      than the maximum power in `df_power`, please
      decrease it!\033[0m"
      )
    )
  }

  # scale down the flow cell capacity according to reads_valid_rate
  flow_capacities <- flow_capacities * reads_valid_rate

  # scale up the cell number according to cells_used_rate
  df_power$N <- df_power$N / cells_used_rate

  # check flow constraints
  if (max(flow_capacities) < min(df_power$N * df_power$R)) {
    warning(
      strwrap(
        "\033[31m The maximum of `flow_capacities` is smaller
      than the minimum of `N*R` in `df_power`, the result
      may be biased!\033[0m"
      )
    )
  }

  # scam model
  model_scam <-
    scam::scam(power ~ s(N, k = 10, bs = 'mpi') + s(R, k = 10, bs = 'mpi'),
               data = df_power)

  df_flow_list <- list()

  for (i in seq_along(flow_capacities)) {
    R <- flow_capacities[i] / unique(df_power$N)
    df_flow_list[[i]] <- data.frame(
      N = unique(df_power$N),
      R = R,
      power = predict(model_scam,
                      newdata = data.frame(
                        N = unique(df_power$N),
                        R = R
                      ))
    ) %>%
      dplyr::filter(N < max(df_power$N) &
                      R < max(df_power$R)) %>%
      dplyr::mutate(flow_cell = i,
                    cost = library_costs + flow_costs[i])
  }

  df_flow <- do.call(rbind, df_flow_list) %>%
    dplyr::mutate(flow_cell = factor(flow_cell))

  temp <- function(y) {
    function(x)
      y / x
  }

  capacties_afforable <-
    (budget - library_costs) / (min(flow_costs / flow_capacities, 1))

  curves_flow_shared <-
    purrr::map(c(flow_capacities, capacties_afforable), ~ temp(.x))

  # curves_flow_ind <- map(flow_capacities, ~ temp(.x))

  converted_capacities <-
    scales::label_number(scale_cut = scales::cut_short_scale())(flow_capacities)
  levels(df_flow$flow_cell) <- converted_capacities

  #------------------------Individual design------------------------#
  ind_available_design <- df_flow %>%
    dplyr::filter(power >= power_threshold,
                  cost <= budget)

  if (nrow(ind_available_design) == 0) {
    warning(
      strwrap(
        "\033[33m Missing optimal indepedent design with current budget
        and power restriction, please increase budget or lower
        power restriction!\033[0m"
      )
    )
  }

  ind_optimal_power <- ind_available_design %>%
    dplyr::arrange(dplyr::desc(power)) %>%
    dplyr::slice_head() %>%
    dplyr::mutate(type = "power")

  ind_optimal_cost <- ind_available_design %>%
    dplyr::arrange(cost, dplyr::desc(power)) %>%
    dplyr::slice_head() %>%
    dplyr::mutate(type = "cost")

  ind_optimal <- rbind(ind_optimal_power, ind_optimal_cost)

  # +2 because of the warning of minimum n is 3
  colors_ind <-
    RColorBrewer::brewer.pal(n = length(flow_capacities) + 3, name = "Set1")[seq_len(1 + length(flow_capacities))]
  names(colors_ind) <- c(converted_capacities, "Budget")

  type_values <- c("power" = 8, "cost" = 9)

  p_design_ind <- ggplot2::ggplot(df_flow,
                                  ggplot2::aes(N, R)) +
    ggplot2::geom_point(shape = 21,
                        size = 3,
                        ggplot2::aes(color = flow_cell, fill = power)) +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "red"),
      breaks = seq(0, 1, length.out = 6),
      limits = c(0, 1),
      values = c(0, 1),
      name = "Power"
    ) +
    purrr::map2(
      curves_flow_shared,
      c(converted_capacities, "Budget"),
      ~ ggplot2::stat_function(fun = .x,
                               ggplot2::aes(color = .y))
    ) +
    ggplot2::scale_color_manual(
      values = colors_ind,
      breaks = c(converted_capacities, "Budget"),
      name = "Flow Cell Capacity"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(suffix = " K", scale = 1e-3),
      limits = c(0, max(df_power$R))
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_number(suffix = " K", scale = 1e-3),
      limits = c(0, max(df_power$N))
    ) +
    ggplot2::labs(x = "Target number of cells",
                  y = "Target reads per cell") +
    ggplot2::geom_point(
      data = ind_optimal[!is.na(ind_optimal$N), ],
      ggplot2::aes(N, R, shape = type),
      size = 4,
      color = "purple"
    ) +
    ggplot2::scale_shape_manual(values = type_values,
                                name = "Optimal Design") +
    ggplot2::theme_bw()

  p_power_cost_ind <- ggplot2::ggplot(df_flow,
                                      ggplot2::aes(cost, power)) +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "red"),
      breaks = seq(0, 1, length.out = 6),
      limits = c(0, 1),
      values = c(0, 1),
      name = "Power"
    ) +
    ggplot2::geom_point(shape = 21,
                        size = 3,
                        ggplot2::aes(color = flow_cell, fill = power)) +
    ggplot2::scale_color_manual(
      values = colors_ind,
      breaks = c(converted_capacities, "Budget"),
      name = "Flow Cell Capacity"
    ) +
    ggplot2::labs(x = "Cost",
                  y = "Power") +
    ggplot2::geom_hline(yintercept = power_threshold, col = "red") +
    ggplot2::geom_vline(xintercept = budget, col = "blue") +
    ggplot2::geom_text(
      x = budget,
      y = min(df_flow$power) + diff(range(df_flow$power)) / 10,
      label = "Budget",
      col = "blue"
    ) +
    ggplot2::geom_text(
      x = max(df_flow$cost, na.rm = TRUE) - diff(range(df_flow$cost) / 5),
      y = power_threshold,
      label = "Threshold",
      col = "red"
    ) +
    ggplot2::geom_point(
      data = ind_optimal[!is.na(ind_optimal$N), ],
      ggplot2::aes(cost, power, shape = type),
      size = 4,
      color = "purple"
    ) +
    ggplot2::scale_shape_manual(values = type_values,
                                name = "Optimal Design") +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw()

  p_ind <- patchwork::wrap_plots(
    p_design_ind +
      p_power_cost_ind +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
  )

  #------------------------Shared design------------------------#
  df_power_shared <- df_power %>%
    dplyr::mutate(cost = library_costs + min(flow_costs / flow_capacities, 1) * N * R)
  shared_flow_cells_index <- which.min(flow_costs / flow_capacities)
  if (shared_flow_cells_index > length(flow_capacities)) {
    share_available_design <- NULL
    share_optimal_cost <- NULL
    share_optimal_power <- NULL
    p_share <- NULL
  } else {
    flow_cell_used <- converted_capacities[shared_flow_cells_index]
    df_power_shared$flow_cell <- flow_cell_used
    share_available_design <- df_power_shared %>%
      dplyr::filter(power >= power_threshold,
                    cost <= budget)

    if (nrow(share_available_design) == 0) {
      warning(
        strwrap(
          "\033[33m Missing optimal share design with current budget
        and power restriction, please increase budget or lower
        power restriction!\033[0m"
        )
      )
    }

    share_optimal_power <- share_available_design %>%
      dplyr::arrange(dplyr::desc(power)) %>%
      dplyr::slice_head() %>%
      dplyr::mutate(type = "power")

    share_optimal_cost <- share_available_design %>%
      dplyr::arrange(cost, dplyr::desc(power)) %>%
      dplyr::slice_head() %>%
      dplyr::mutate(type = "cost")

    share_optimal <- rbind(share_optimal_power,
                           share_optimal_cost)

    p_design_share <- ggplot2::ggplot(df_power_shared,
                                      ggplot2::aes(N, R)) +
      ggplot2::geom_point(shape = 21,
                          size = 3,
                          ggplot2::aes(fill = power)) +
      ggplot2::scale_fill_gradientn(
        colors = c("white", "red"),
        breaks = seq(0, 1, length.out = 6),
        limits = c(0, 1),
        values = c(0, 1),
        name = "Power"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(suffix = " K", scale = 1e-3),
        limits = c(0, max(df_power$R))
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::label_number(suffix = " K", scale = 1e-3),
        limits = c(0, max(df_power$N))
      ) +
      purrr::map2(
        curves_flow_shared,
        c(converted_capacities, "Budget"),
        ~ ggplot2::stat_function(fun = .x,
                                 ggplot2::aes(color = .y))
      ) +
      ggplot2::scale_color_manual(
        values = colors_ind,
        breaks = c(converted_capacities, "Budget"),
        name = "Flow Cell Capacity"
      ) +
      ggplot2::geom_point(
        data = share_optimal,
        ggplot2::aes(N, R, shape = type),
        size = 4,
        color = "purple"
      ) +
      ggplot2::scale_shape_manual(values = type_values,
                                  name = "Optimal Design") +
      ggplot2::labs(x = "Target number of cells",
                    y = "Target reads per cell",
                    shape = "Optimal Design") +
      ggplot2::theme_bw()

    p_power_cost_share <- ggplot2::ggplot(df_power_shared,
                                          ggplot2::aes(cost, power, fill = power)) +
      ggplot2::scale_fill_gradientn(
        colors = c("white", "red"),
        breaks = seq(0, 1, length.out = 6),
        limits = c(0, 1),
        values = c(0, 1),
        name = "Power"
      ) +
      ggplot2::geom_point(shape = 21,
                          size = 3,
                          ggplot2::aes(color = flow_cell)) +
      ggplot2::labs(x = "Cost",
                    y = "Power") +
      ggplot2::scale_color_manual(
        values = colors_ind,
        breaks = c(converted_capacities, "Budget"),
        name = "Flow Cell Capacity"
      ) +
      ggplot2::geom_hline(yintercept = power_threshold, col = "red") +
      ggplot2::geom_vline(xintercept = budget, col = "blue") +
      ggplot2::geom_text(
        x = budget,
        y = min(df_power_shared$power) + diff(range(df_power_shared$power)) / 10,
        label = "Budget",
        col = "blue"
      ) +
      ggplot2::geom_text(
        x = max(df_power_shared$cost, na.rm = TRUE) - diff(range(df_power_shared$cost) / 5),
        y = power_threshold,
        label = "Power Threshold",
        col = "red"
      ) +
      ggplot2::geom_point(
        data = share_optimal,
        ggplot2::aes(cost, power, shape = type),
        size = 4,
        color = "purple"
      ) +
      ggplot2::scale_shape_manual(values = type_values,
                                  name = "Optimal Design") +
      ggplot2::theme_bw() +
      ggplot2::guides(color = "none")

    p_share <- patchwork::wrap_plots(
      p_design_share +
        p_power_cost_share +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    )
  }

  #------------------------Return------------------------#
  output <- list(
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
