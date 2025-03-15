#' @export
rsizes <- function(n_items, n_clusters, alpha, beta) {
  .Call(.rsizes, n_items, n_clusters, alpha, beta)
}

#' @export
dsizes <- function(x, alpha, beta, log = FALSE) {
  .Call(.dsizes, x, alpha, beta, log)
}

#' @export
rpartition <- function(n_items, n_clusters, alpha, beta) {
  .Call(.rpartition, n_items, n_clusters, alpha, beta)
}

#' @export
cfsc <- function(n_items, n_clusters, w) {
  .Call(.count_for_size_configuration, n_items, n_clusters, w)
}
