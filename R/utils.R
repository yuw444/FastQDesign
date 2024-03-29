#' Jaccard Index
#'
#' Calculate the Jaccard index between two sets
#' @param a an vector or list
#' @param b an vector or list
#' @return scalar, jaccard index
#' @examples
#' set.seed(926)
#' a <- sample(1:100, 50)
#' b <- sample(1:100, 50)
#' JaccardIndex(a,b)
#' @export
JaccardIndex <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(round(intersection / union, 3))
}

#' Positive Predictive Value(PPV)
#'
#' Calculate the PPV index between two sets
#' @param a an vector or list, predicted value
#' @param b an vector or list, true value
#' @return scalar, PPV index
#' @examples
#' set.seed(926)
#' a <- sample(1:100, 50)
#' b <- sample(1:100, 50)
#' PPV(a,b)
#' @export
#'
PPV <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  intersection <- length(intersect(a, b))
  return(round(intersection / length(a), 5))
}

#' Entropy
#'
#' Calculate the entropy from the probability mass function(PMF)
Entropy <- function(x) {
  out <- 0
  for (e in x) {
    if (e != 0) {
      out <- out - e * log2(e)
    }
  }
  return(out)
}

#' Conditional Shannon Entropy
#'
#' Calculate Conditional Shannon Entropy, weighted by its proportion
#' on the condition
#'
#' @param condition an vector, the condition
#' @param feature an vector, the interested feature
#' @return a scalar
#'
#' @examples
#' set.seed(926)
#' condition <- sample(1:100, 50)
#' feature <- sample(1:100, 50)
#' ConditionalShannonEntropy(condition,feature)
#' @export
#'
ConditionalShannonEntropy <- function(condition, feature) {
  temp <- table(feature, condition)

  temp <- matrix(temp, ncol = ncol(temp), dimnames = dimnames(temp))

  tempRowSums <- apply(temp, 1, sum)

  tt <- temp / tempRowSums

  rowProp <- tempRowSums / sum(temp)

  out <- sum(apply(tt, 1, Entropy) * rowProp)

  return(out)
}


CalculatePowerIndexes <- function(ListReference,
                                  ListDownsample,
                                  ...) {
  require(Seurat)
  require(dplyr)
  require(tibble)
  require(mclust)
  require(igraph)

  cluster_info <- data.frame(
    Id = Cells(ListReference$Seurat),
    Cluster = ListReference$Seurat$seurat_clusters
  )

  cluster_info_downsample <- data.frame(
    Id = Cells(ListDownsample$Seurat),
    Cluster = ListDownsample$Seurat$seurat_clusters
  )

  matched_cluster_info <- cluster_info %>%
    left_join(., cluster_info_downsample, by = c("Id" = "Id")) %>%
    select(Cluster.y) %>%
    unlist()

  geneList_cluster_reference <- MarkerGeneFilter(ListReference$MarkersByCluster, ...)
  geneList_condition_reference <- MarkerGeneFilter(ListReference$MarkersByCondition, ...)
  geneList_cluster_downsample <- MarkerGeneFilter(ListDownsample$MarkersByCluster, ...)
  geneList_condition_downsample <- MarkerGeneFilter(ListDownsample$MarkersByCondition, ...)

  ARI <- adjustedRandIndex(cluster_info$Cluster, matched_cluster_info)
  ConditionalEntropy <- ConditionalShannonEntropy(cluster_info$Cluster, matched_cluster_info)
  JaccardIndex_cluster <- JaccardIndex(geneList_cluster_reference, geneList_cluster_downsample)
  JaccardIndex_condition <- JaccardIndex(geneList_condition_reference, geneList_condition_downsample)

  return(list(
    ARI = ARI,
    ConditionalEntropy = ConditionalEntropy,
    JaccardIndex_cluster = JaccardIndex_cluster,
    JaccardIndex_condition = JaccardIndex_condition
  ))
}

PowerCalPipeline <- function(downsample,
                             ListReference,
                             fixedNCluster = TRUE,
                             resolution,
                             n.clusters,
                             verbose = FALSE,
                             ...) {
  require(Seurat)
  require(glmGamPoi)
  downsample <- SCTransform(
    object = downsample,
    method = "glmGamPoi",
    verbose = verbose
  )
  downsample <- RunPCA(downsample, verbose = verbose)
  downsample <- FindNeighbors(downsample, verbose = verbose)

  if (fixedNCluster == TRUE) {
    downsample <- FixedNumClusters(downsample, n.clusters)
  } else {
    downsample <- FindClusters(downsample, resolution = resolution)
  }

  downsample <- RunUMAP(downsample, dims = 1:30, verbose = verbose)

  markers_downsample <- FindAllMarkers(
    downsample,
    min.pct = 0.1,
    logfc.threshold = 0,
    min.cells.feature = 0,
    min.cells.group = 0,
    return.thresh = 0.1,
    slot = "data",
    verbose = verbose
  )

  markers_downsample_condition <- FindAllMarkersByCondition(
    downsample,
    condition = "orig.ident",
    min.pct = 0.1,
    logfc.threshold = 0,
    min.cells.feature = 0,
    min.cells.group = 0,
    return.thresh = 0.1,
    slot = "data",
    verbose = verbose
  )

  downsample_list <- list(
    Seurat = downsample,
    MarkersByCluster = markers_downsample,
    MarkersByCondition = markers_downsample_condition
  )

  downsample_list$Indexes <- CalculatePowerIndexes(
    ListReference,
    downsample_list, ...
  )

  return(downsample_list)
}

ExtractPowerIndex <- function(file) {
  temp <- readRDS(file = file)
  return(unlist(temp$Indexes))
}
