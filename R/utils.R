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
#'
#'
ExtractMTX <- function(file_path){

  metadata_lines <- readLines(file_path, n = 20)
  metadata_lines[5] <- paste0(metadata_lines[5], ",", sep="")
  metadata_lines[6] <- "\t\"parent_bam\": \"/scratch/g/chlin/Yu/AIBM/sub/AI/outs/possorted_genome_bam.bam\","

  # Extract metadata lines containing '%metadata_json'
  metadata_json_lines <- gsub("^%", "", metadata_json_lines)
  metadata_json <- paste(metadata_json, collapse = "")
  metadata_json <- regmatches(metadata_json, regexpr("\\{.*\\}", metadata_json))

  # Convert JSON string to a list
  metadata_list <- jsonlite::fromJSON(metadata_json)

  return(metadata_list)
}
