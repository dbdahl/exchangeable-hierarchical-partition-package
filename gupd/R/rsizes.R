#' @export
rsizes <- function(n_items, n_clusters, alpha, beta) {
  .Call(.rsizes, n_items, n_clusters, alpha, beta)
}
