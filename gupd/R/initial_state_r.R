#' @export
new_from_vec <- function(x) {
  .Call(.new_from_vec, x)
}

#' @export
new_max_entropy <- function(n_items, n_clusters) {
  .Call(.new_max_entropy, n_items, n_clusters)
}

#' @export
size_configuration_to_r <- function(x) {
  .Call(.size_configuration_to_r, x)
}

#' @export
size_configuration_available <- function(x, index) {
  .Call(.size_configuration_available, x, index)
}

#' @export
size_configuration_redistribute <- function(x, index, n) {
  .Call(.size_configuration_redistribute, x, index, n)
}
