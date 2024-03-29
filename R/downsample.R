#' Down Sample a Seurat object
#'
#' Down sample a seurat object with given `rate_cells` and `rate_umis`
#'
#' @import Seurat
#' @param seu A Seurat object
#' @param rate_cells The proportion of cell to sample
#' @param rate_umis The proportion of UMIs to sample
#' @param seed A seed for random number generation
#' @param enable_PCR Whether to consider the duplication of UMIs,
#'                   default is `TRUE`, use negative binomial distribution
#'                   to simulate the number of copies of each UMI
#' @param nb_size If `enable_PCR=TRUE`, it is used in \link{rnbinom}
#' @param nb_prob If `enable_PCR=TRUE`, it is used in \link{rnbinom}
#'
#' @export

DownSample <- function(
    seu,
    rate_cells = 0.3,
    rate_umis = 0.3,
    seeds = 926,
    enable_PCR = FALSE,
    nb_size = 2,
    nb_prob = 0.2) {
  set.seed(seeds * 1000 + rate_cells * 100 + rate_umis * 10)

  n_obs <- ncol(seu)
  if (rate_umis == 1) {
    j <- sort(sample(1:n_obs, ceiling(n_obs * rate_cells)))
    seu_down <- subset(seu, cells = j)
    Seurat::DefaultAssay(seu_down) <- "RNA"
  } else if (rate_cells == 1) {
    umis <- seu@assays[["RNA"]]@counts@x
    n_umis <- length(umis)
    if (enable_PCR) {
      cat("Enabling copy number variation of each UMI...\nDownsampling Starting...\n")
      pb <- txtProgressBar(min = 1, max = n_umis, style = 3)
      umis_new <- sapply(1:n_umis, function(x) {
        umis_copy_number <- rnbinom(umis[x], nb_size, nb_prob)
        setTxtProgressBar(pb, x)
        return(sum(rbinom(umis[x], umis_copy_number, rate_umis) != 0))
      })
      cat("\n")
    } else {
      umis_new <- rbinom(n_umis, umis, rate_umis)
    }
    seu@assays[["RNA"]]@counts@x <- as.numeric(umis_new)
    seu_down <- seu
    seu_down@assays[["RNA"]]@counts <- as(as.matrix(seu@assays[["RNA"]]@counts), "dgCMatrix")
    seu_down$nCount_RNA <- colSums(seu_down@assays[["RNA"]]@counts)
    seu_down$nFeature_RNA <- diff(seu_down@assays[["RNA"]]@counts@p)
    Seurat::DefaultAssay(seu_down) <- "RNA"
  } else {
    j <- sort(sample(1:n_obs, ceiling(n_obs * rate_cells)))
    temp <- subset(seu, cells = j)
    umis <- temp@assays[["RNA"]]@counts@x
    n_umis <- length(umis)
    if (enable_PCR) {
      cat("Enabling copy number variation of each UMI...\nDownsampling Starting...\n")
      pb <- txtProgressBar(min = 1, max = n_umis, style = 3)
      umis_new <- sapply(1:n_umis, function(x) {
        umis_copy_number <- rnbinom(umis[x], nb_size, nb_prob)
        setTxtProgressBar(pb, x)
        return(sum(rbinom(umis[x], umis_copy_number, rate_umis) != 0))
      })
      cat("\n")
    } else {
      umis_new <- rbinom(n_umis, umis, rate_umis)
    }
    temp@assays[["RNA"]]@counts@x <- as.numeric(umis_new)
    seu_down <- temp
    seu_down@assays[["RNA"]]@counts <- as(as.matrix(temp@assays[["RNA"]]@counts), "dgCMatrix")
    seu_down$nCount_RNA <- colSums(seu_down@assays[["RNA"]]@counts)
    seu_down$nFeature_RNA <- diff(seu_down@assays[["RNA"]]@counts@p)
    Seurat::DefaultAssay(seu_down) <- "RNA"
  }

  return(seu_down)
}


