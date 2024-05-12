#' Down Sample a Seurat object
#'
#' Down sample a seurat object with given `rate_cells` and `rate_umis`
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Matrix colSums
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

DownSample <- function(seu,
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
      pb <- txtProgressBar(min = 1,
                           max = n_umis,
                           style = 3)
      umis_new <- sapply(1:n_umis, function(x) {
        umis_copy_number <- rnbinom(umis[x], nb_size, nb_prob)
        setTxtProgressBar(pb, x)
        return(sum(rbinom(
          umis[x], umis_copy_number, rate_umis
        ) != 0))
      })
      cat("\n")
    } else {
      umis_new <- rbinom(n_umis, umis, rate_umis)
    }
    seu@assays[["RNA"]]@counts@x <- as.numeric(umis_new)
    seu_down <- seu
    seu_down@assays[["RNA"]]@counts <-
      as(as.matrix(seu@assays[["RNA"]]@counts), "dgCMatrix")
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
      pb <- txtProgressBar(min = 1,
                           max = n_umis,
                           style = 3)
      umis_new <- sapply(1:n_umis, function(x) {
        umis_copy_number <- rnbinom(umis[x], nb_size, nb_prob)
        setTxtProgressBar(pb, x)
        return(sum(rbinom(
          umis[x], umis_copy_number, rate_umis
        ) != 0))
      })
      cat("\n")
    } else {
      umis_new <- rbinom(n_umis, umis, rate_umis)
    }
    temp@assays[["RNA"]]@counts@x <- as.numeric(umis_new)
    seu_down <- temp
    seu_down@assays[["RNA"]]@counts <-
      as(as.matrix(temp@assays[["RNA"]]@counts), "dgCMatrix")
    seu_down$nCount_RNA <- Matrix::colSums(seu_down@assays[["RNA"]]@counts)
    seu_down$nFeature_RNA <- diff(seu_down@assays[["RNA"]]@counts@p)
    Seurat::DefaultAssay(seu_down) <- "RNA"
  }

  return(seu_down)
}


#' Prepare a sample
#'
#' Prepare a sample for the comparison
#'
#' @importFrom Seurat SCTransform RunPCA RunUMAP FindNeighbors FindAllMarkers FindClusters
#' @import glmGamPoi
#' @importFrom monocle3 pseudotime
#' @importFrom SummarizedExperiment colData
#' @param seu A \code{Seurat} object
#' @param cluster_id A A character, the condition colname in \code{Seurat@meta.data}
#' @param n_clusters A scalar, the desired number of clusters
#' @param use_default_res Whether to use the default resolution in \code{Seurat::FindClusters}
#' @param condition A character, the condition colname in \code{Seurat@meta.data}
#' @param cell_3d_embedding Calculate 3d embedding for the input \code{seu}
#' @param pca_used Used in \code{Seurat::FindNeighbors} and \code{Seurat::RunUMAP}
#' @param vars_to_regress A character string, used in \code{Seurat::SCTransform}
#' @param use_partition Whether to \code{use_partition} when \code{monocle3::learn_graph}
#' @param close_loop Logic, a parameter in \code{monocle3::learn_graph}
#' @param learn_graph_control List, a parameter in \code{monocle3::learn_graph}
#' @param interactive Logic, whether to select root cell interactively when considering cell fate
#' @param root_cells_ref An vector of cell ids can be used as the root of cell trajectory
#' @return list Seurat with extra meta.data features,
#'  data.frames of DE genes by cluster and condition
#'
#' @export
#'
SamplePrep <- function(seu,
                       cluster_id = NA,
                       n_clusters = NA,
                       use_default_res = TRUE,
                       condition = NA,
                       cell_3d_embedding = FALSE,
                       pca_used = 1:30,
                       vars_to_regress = NULL,
                       use_partition = FALSE,
                       close_loop = FALSE,
                       learn_graph_control = NULL,
                       interactive = FALSE,
                       root_cells_ref = NA,
                       verbose = FALSE,
                       ...) {

  if(!is.na(cluster_id) & (!is.na(n_clusters) | use_default_res))
  {
    warning("only applied `cluster_id`, and ignored `n_clusters` and `use_default_res.`\n")
    n_clusters <- NA
    use_default_res <- FALSE
  }

  if (!is.na(n_clusters) | use_default_res)
  {
    cat("\n###When n_cluster is specified, SCTransform is applied automatically!\n")
    seu <- Seurat::SCTransform(object = seu,
                               method = "glmGamPoi",
                               vars.to.regress = vars_to_regress,
                               verbose = verbose)
    seu <- Seurat::RunPCA(seu, verbose = verbose)
    seu <- Seurat::FindNeighbors(seu, dims = pca_used, verbose = verbose)
    if (!is.na(n_clusters))
    {
      if(use_default_res){
        warning("Ignored `use_default_res.`\n")
      }
      cat("\n###Finding the desired number of clusters ...\n")
      seu <- FixedNumClusters(seu, n_clusters)
    } else {
      seu <- Seurat::FindClusters(seu)
    }
    seu <- Seurat::RunUMAP(seu, dims = pca_used, verbose = verbose)
  }

  if (!is.na(root_cells_ref[1]) | interactive) {
    cat("\n### Finding the pseudotime ....\n")
    temp_cds <- FastQDesign::FindPseudotime(
      seu,
      redo_sctransform = FALSE,
      vars_to_regress = vars_to_regress,
      use_partition = use_partition,
      close_loop = close_loop,
      learn_graph_control = learn_graph_control,
      interactive = interactive,
      root_cells_ref = root_cells_ref
    )
    seu$pseudotime <- monocle3::pseudotime(temp_cds)
    seu$root_cells <- SummarizedExperiment::colData(temp_cds)$root_cells
  }

  cat("\n### Finding all cluster markers ....\n")
  if(!is.na(cluster_id))
    Idents(seu) <- cluster_id

  df_marker_cluster <- Seurat::FindAllMarkers(seu, ...)

  df_marker_condition <- NULL

  if (!is.na(condition)) {
    cat("\n### Finding all condition markers ....\n")
    df_marker_condition <- FindAllMarkersByCondition(seu,
                                                     condition = condition,
                                                     ...)
  }

  seu$pseudotime <- NA

  seu$umap_1 <- seu$umap_2 <- seu$umap_3 <- NA
  if (cell_3d_embedding) {
    cat("\n### Finding the cells 3d embeddings ....\n")
    suppressWarnings({
      temp <-
        Seurat::RunUMAP(seu,
                        dims = pca_used,
                        n.components = 3,
                        verbose = FALSE)
      umap3 <- Seurat::Embeddings(temp, "umap")
    })
    seu$umap_1 <- umap3[, 1]
    seu$umap_2 <- umap3[, 2]
    seu$umap_3 <- umap3[, 3]
  }

  downsample_list <- list(
    Seurat = seu,
    MarkersByCluster = df_marker_cluster,
    MarkersByCondition = df_marker_condition
  )

  return(downsample_list)
}

#' Match Downsample to Reference
#'
#' Match the features of downsample to the reference, cluster membership,
#' differential expressed genes by cluster and condition, pseudotime,
#' cell lower dimension embedding
#'
#' @importFrom magrittr %>%
#' @param reference_list The list structure same with \code{SamplePrep} returns
#' @param downsample_list The list structure same with \code{SamplePrep} returns
#' @param ... Parameters that pass to \code{MarkerGeneFilter}
#'
#' @return List of matched features
#' @export
SampleMatch <- function(reference_list,
                        downsample_list,
                        ...) {
  df_cluster_ref <- data.frame(Id = Cells(reference_list$Seurat),
                               reference_list$Seurat@meta.data[, c("seurat_clusters",
                                                                   "pseudotime",
                                                                   "umap_1",
                                                                   "umap_2",
                                                                   "umap_3")])

  df_cluster_ds <- data.frame(Id = Cells(downsample_list$Seurat),
                              downsample_list$Seurat@meta.data[, c("seurat_clusters",
                                                                   "pseudotime",
                                                                   "umap_1",
                                                                   "umap_2",
                                                                   "umap_3")])

  df_cluster_match <- df_cluster_ds %>%
    dplyr::left_join(., df_cluster_ref, by = c("Id" = "Id"))

  genes_cluster_ref <-
    MarkerGeneFilter(reference_list$MarkersByCluster, ...)
  genes_cluster_ds <-
    MarkerGeneFilter(downsample_list$MarkersByCluster, ...)

  genes_condition_ref <- genes_condition_ds <- NA
  if (!is.null(reference_list$MarkersByCondition) &
      !is.null(downsample_list$MarkersByCondition))
  {
    genes_condition_ref <-
      MarkerGeneFilter(reference_list$MarkersByCondition, ...)
    genes_condition_ds <-
      MarkerGeneFilter(downsample_list$MarkersByCondition, ...)
  }

  return(
    list(
      cluster_match = df_cluster_match,
      genes_cluster = list(ref = genes_cluster_ref,
                           ds = genes_cluster_ds),
      genes_condition = list(ref = genes_condition_ref,
                             ds = genes_condition_ds)

    )
  )
}


