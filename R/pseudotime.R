#' Root Node Selection for \code{cell_data_set}
#'
#' Select the root node for \code{cell_data_set} based on the given \code{root_cells_ref}
#' @import monocle3 igraph
#' @param cds A \code{cell_data_set} object from \code{monocle3}
#' @param root_cells_ref An vector of cell ids can be used as the root of cell trajectory
#' @return A character, rootnode id
#' @export
RootNodeSelect <- function(cds,
                           root_cells_ref = NA) {
  cell_closet_vertex <-
    cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]][, 1]
  leaf_node <-
    which(igraph::degree(monocle3::principal_graph(cds)[["UMAP"]]) == 1)

  if (any(is.na(root_cells_ref))) {
    # If root cells are not provided, choose the biggest node as root node
    candidate_root_1 <-
      names(which.max(table(cell_closet_vertex[cell_closet_vertex %in% leaf_node])))
    # In case, tie happens, then select the first one
    return(candidate_root_1[1])
  }

  candidate_root_1 <-
    names(which.max(table(cell_closet_vertex[(cell_closet_vertex %in% leaf_node) &
                                               (names(cell_closet_vertex) %in% root_cells_ref)])))

  # if there is no root cell in the leaf nodes, then select the leaf node with the most cells
  if (!is.null(candidate_root_1)) {
    return(candidate_root_1[1])
  }

  candidate_root_2 <-
    names(which.max(table(cell_closet_vertex[cell_closet_vertex %in% leaf_node])))
  warning(
    "None of provided 'root_cells_ref' is in the leaf node, choose the biggest leaf node as the root node!"
  )
  return(candidate_root_2[1])
}

#' Find Pseudotime for \code{Seurat} Object
#'
#' Find the pseudotime of a \code{Seurat} object based on its \code{UMAP}.
#'
#' @import monocle3 Seurat
#' @param seu A \code{Seurat} object
#' @param root_cells_ref An vector of cell ids can be used as the root of cell trajectory
#' @param min_branch_len The minimum branch length in the cell trajectory, the same parameter in \code{monocle3::learn_graph}
#' @param redo_sctransform Whether to redo the \code{Seurat::SCTransform}
#' @return A \code{cell_data_set} object with pseudotime as feature
#' @export
#'
FindPseudotime <- function(seu,
                           root_cells_ref = NA,
                           min_branch_len = 10,
                           redo_sctransform = FALSE) {
  require(monocle3)
  require(SeuratWrappers)
  require(tibble)

  if (redo_sctransform) {
    seu <- Seurat::SCTransform(seu, verbose = TRUE)
    seu <- Seurat::RunPCA(seu, verbose = TRUE)
    seu <-
      Seurat::RunUMAP(seu,
                      reduction = "pca",
                      dims = 1:30,
                      verbose = FALSE)
  }
  cds <- Seurat::as.CellDataSet(seu)

  cds <-
    monocle3::cluster_cells(cds = cds, reduction_method = "UMAP")
  cds <- monocle3::learn_graph(
    cds,
    use_partition = FALSE,
    learn_graph_control = list(minimal_branch_len = min_branch_lens)
  )

  cell_closet_vertex <-
    cds@principal_graph_aux@listData[["UMAP"]][["pr_graph_cell_proj_closest_vertex"]]

  root_nodes <- RootNodeSelect(cds, root_cells_ref)

  cds <- monocle3::order_cells(cds,
                               root_pr_nodes = paste0("Y_", root_nodes))

  cds$cell_closet_vertex <- cell_closet_vertex
  cds$root_cells_ref <- as.factor(cell_closet_vertex == root_nodes)
  cds$pseudotime <-
    cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

  return(cds)
}
