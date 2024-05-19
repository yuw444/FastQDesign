#' Customize Cluster Number
#'
#' Customize the number of cluster in the Seurat object
#'
#' @param seu A \code{Seurat} object
#' @param n_clusters A scalar, the desired number of clusters
#' @param verbose Default is \code{\link{FALSE}}
#'
#' @export

FixedNumClusters <- function(seu,
                             n_clusters,
                             verbose = FALSE) {
  res_ <- uniroot(
    function(res) {
      nlevels(FindClusters(seu, resolution = res, verbose = verbose)$seurat_clusters) - n_clusters
    },
    interval = c(0, 4), extendInt = "yes"
  )$root

  seu <- FindClusters(seu, resolution = res_, verbose = verbose)

  return(seu)
}

#' Find the optimal number of cell clusters for Seurat object
#'
#' Find the optimal number of cell clusters using \code{fpc::prediction.strenth}
#'
#' @importFrom fpc prediction.strength
#'
#'
#' @param seu A seurat objet
#' @param reduction Reduction to use as input for \code{fpc:prediction.strength}
#'
#' @return An object of class \code{predstr}, please see \code{fpc} for the detail
#' @export

OptimalNumClusters <- function(
    seu,
    reduction = c("pca", "umap")
){

  assertthat::assert_that(any(class(seu) == "Seurat"), msg = "Please make sure seu is an Seurat object")
  cell_embeddings <- Seurat::Embeddings(seu, reduction = reduction)
  if(reduction == "pca")
  {
    save_commands <- Seurat::Command(seu)
    which_pca <- which(sapply(save_commands, function(x) grepl(".pca", x)))[1]
    dims <- Seurat::Command(seu, command = save_commands[which_pca], "dims")
    cell_embeddings <- cell_embeddings[, dims]
  }

  return(fpc::prediction.strength(cell_embeddings))

}
