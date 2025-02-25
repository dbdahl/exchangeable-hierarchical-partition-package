#' Generalized Uniform Partition Distribution
#'
#' A function implementing the probability mass function is provided based on
#' the supplied probabilty vector and number of items.
#'
#' @param prob An vector of probability weights, where element \eqn{i} is the
#' unnormalized marginal probabilty of partitions with \eqn{i} clusters.  If
#' the length of this vector is less than the number of items, zero is implied
#' for all remaining elements.
#' @param n_items An integer for the number of items in the partition.
#'
#' @return A function to evaluate the Generalized Uniform Partition Distribution.
#' @export
#' @examples
#' f <- make_gupd(c(1, 10, 5), 6)
#' f(c(1, 1, 2, 2, 3, 3))
#' f(3)
#'
make_gupd <- function(prob, n_items) {
  x <- .Call(.make_gupd, prob, n_items)
  function(k) {
    msg <- function() stop(sprintf("'k' should be: 1. an integer greater than 0, or 2: a vector of cluster labels of length %s", n_items))
    if (length(k) == 0) msg()
    if (length(k) == n_items) k <- length(unique(k))
    if (length(k) != 1) msg()
    if (k <= length(x)) x[k] else 0.0
  }
}
