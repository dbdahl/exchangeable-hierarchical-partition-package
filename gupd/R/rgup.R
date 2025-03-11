#' @export
rgup <- function(n_items, n_clusters, alpha, beta) {
  .Call(.rgup, n_items, n_clusters, alpha, beta)
}
