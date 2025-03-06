#' Generalized Uniform Partition Distribution
#'
#' A function implementing the probability mass function is provided based on
#' the supplied probabilty vector and number of items.
#'
#' @param n_items An integer for the number of items in the partition.
#' @param n_clusters An integer for the number of clusters in the partition.
#' @param n_samples An integer for the number of sampled partitions.
#' @param a Parameter
#'
#' @return A function to evaluate the Generalized Uniform Partition Distribution.
#' @export
#' @examples
#' experiment(10, 3, 5, 1)
#'
experiment <- function(n_items, n_clusters, n_samples, concentration) {
  .Call(.experiment, n_items, n_clusters, n_samples, concentration)
}
