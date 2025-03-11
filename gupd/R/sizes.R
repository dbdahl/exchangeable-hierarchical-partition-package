#' @export
rsizes <- function(n_items, n_clusters, alpha, beta) {
  .Call(.rsizes, n_items, n_clusters, alpha, beta)
}

#' @export
dsizes <- function(x, alpha, beta, log = FALSE) {
  .Call(.dsizes, x, alpha, beta, log)
}



