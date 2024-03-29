#' Customize Cluster Number
#'
#' Customize the number of cluster in the Seurat object
#'
#' @param seu a Seurat object
#' @param n_clusters scalar, the desired number of clusters
#' @param verbose default is \code{\link{FALSE}}
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

