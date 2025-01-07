#' Find Markers By Condition
#'
#' Find Differential Expressed(DE) markers by Condition
#' @importFrom magrittr %>%
#' @importFrom Seurat FindMarkers Idents
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @param seu A Seurat object
#' @param condition A character name in the `meta.data` of `seu`
#' @param ... The argument could pass to \code{\link{FindMarkers}}
#'
#' @return A data frame includes DE markers
#' @export
FindAllMarkersByCondition <- function(
    seu,
    condition,
    ...) {
  markers_by_condition <- NULL
  k <- 1

  temp <- as.matrix(table(seu$seurat_clusters, seu@meta.data[, condition]))
  temp_qualifer <- apply(temp, 1, function(x) {
    x[1] > 3 & x[2] > 3 & sum(x) > 10
  })

  # Find makers based on Idents
  for (iter in levels(Idents(seu))) {
    if (temp_qualifer[iter]) {
      markers_by_condition[[k]] <- Seurat::FindMarkers(seu,
        ident.1 = seu@meta.data[, condition][1],
        group.by = condition,
        subset.ident = iter,
        ...
      ) %>%
        dplyr::mutate(cluster = iter) %>%
        tibble::rownames_to_column(var = "gene")
    }
    k <- k + 1
  }

  markers_by_condition <- do.call("rbind", markers_by_condition)
  return(markers_by_condition)
}

#' DE Markers Filters
#'
#' Filter the maker genes with series of filters
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter group_by arrange pull
#' @param df_MarkerGene a `data.fram` generated from `FindMarkers`
#' @param pct_1 A minimum cutoff for `pct.1`
#' @param pct_2 A minimum cutoff for `pct.2`
#' @param p_val_adj_ A maximum cutoff for `p_val_adj`
#' @param avg_log2FC_abs A minimum cutoff for `abs(avg_log2FC)`
#'
#' @return An vector of genes that pass through the filters
#' @export
#'
MarkerGeneFilter <- function(df_MarkerGene,
                             pct_1 = 0.2,
                             pct_2 = 0.2,
                             p_val_adj_ = 0.05,
                             avg_log2FC_abs = 0) {
  if (nrow(df_MarkerGene) == 0) {
    return(NA)
  } else {
    geneList <- df_MarkerGene %>%
      dplyr::filter(`p_val_adj` < p_val_adj_ &
        abs(`avg_log2FC`) > avg_log2FC_abs &
        (`pct.1` > pct_1 | `pct.2` > pct_2)) %>%
      dplyr::group_by(`gene`) %>%
      dplyr::arrange(`p_val_adj`) %>%
      dplyr::slice_head() %>%
      dplyr::ungroup() %>%
      dplyr::pull(`gene`)

    return(geneList)
  }
}
