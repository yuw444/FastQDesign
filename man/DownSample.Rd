% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downsample.R
\name{DownSample}
\alias{DownSample}
\title{Down Sample a Seurat object}
\usage{
DownSample(
  seu,
  rate_cells = 0.3,
  rate_umis = 0.3,
  seed = 926,
  enable_PCR = FALSE,
  nb_size = 2,
  nb_prob = 0.2
)
}
\arguments{
\item{seu}{A Seurat object}

\item{rate_cells}{The proportion of cell to sample}

\item{rate_umis}{The proportion of UMIs to sample}

\item{seed}{A seed for random number generation}

\item{enable_PCR}{Whether to consider the duplication of UMIs,
default is `TRUE`, use negative binomial distribution
to simulate the number of copies of each UMI}

\item{nb_size}{If `enable_PCR=TRUE`, it is used in \link{rnbinom}}

\item{nb_prob}{If `enable_PCR=TRUE`, it is used in \link{rnbinom}}
}
\description{
Down sample a seurat object with given `rate_cells` and `rate_umis`
}
