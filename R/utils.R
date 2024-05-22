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


#' Extract metadata from mtx.gz file
#'
#' @importFrom jsonlite fromJSON
#' @param file_path File path to \code{matrix.mtx.gz}
#' @return A list of metadata
#'
#' @export
ExtractMTX <- function (file_path)
{
  con <- file(file_path, "r")
  metadata_lines <- readLines(con, n = 20)
  close(con)
  metadata_json <- gsub("^%", "", metadata_lines)
  metadata_json <- paste(metadata_json, collapse = "")
  metadata_json <- regmatches(metadata_json, regexpr("\\{.*\\}",
                                                     metadata_json))
  metadata_list <- jsonlite::fromJSON(metadata_json)
  c_num <- as.integer(strsplit(metadata_lines[which(grepl("}", metadata_lines)) + 1], " ")[[1]])
  metadata_list$sampled_valid_n_genes <- c_num[1]
  metadata_list$sampled_valid_n_cells <- c_num[2]
  metadata_list$sampled_valid_n_Features <- c_num[3]
  return(metadata_list)
}
