#' Clustering Entropy
#'
#' The clustering entropy is computed.  Let \eqn{p_i} be the relative size of
#' cluster \eqn{i}. The clustering entropy is the negative of the sum of
#' \eqn{p_i * f(p_i)} for all \eqn{i}, where \eqn{f} is the natural logarithm.
#'
#' @param x A vector of cluster labels
#'
#' @return A scalar giving the clustering entropy.
#' @export
#' @examples
#' entropy(c(1, 1, 1, 1, 1))
#' entropy(c(1, 1, 2, 2, 2))
#' entropy(c(1, 2, 3, 4, 5))
#'
entropy <- function(x) {
  .Call(.entropy, x)
}

#' @export
entropy_from_partition <- function(x) {
  .Call(.entropy_from_partition, x)
}

#' @export
entropy_from_cluster_sizes <- function(x) {
  .Call(.entropy_from_cluster_sizes, x)
}
