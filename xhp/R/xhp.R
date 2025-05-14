#' @export
new <- function(n_items, n_clusters_log_weights, cluster_sizes_distribution) {
  .Call(.xhp__new, n_items, n_clusters_log_weights, cluster_sizes_distribution)
}

#' @export
print.xhp <- function(xhp) {
  cat("Pointer to a Generalized Hierarchical Uniform Partition Distribution (xhp)\n")
}

#' @export
sample_cluster_sizes_given_n_clusters <- function(xhp, n_clusters) {
  .Call(.xhp__sample_cluster_sizes_given_n_clusters, xhp, n_clusters)
}

#' @export
sample_n_clusters <- function(xhp) {
  .Call(.xhp__sample_n_clusters, xhp)
}

#' @export
sample_partition <- function(xhp) {
  .Call(.xhp__sample_partition, xhp)
}

#' @export
sample_partition_given_cluster_sizes <- function(xhp, cluster_sizes) {
  .Call(.xhp__sample_partition_given_cluster_sizes, xhp, cluster_sizes)
}

#' @export
sample_partition_given_n_clusters <- function(xhp, n_clusters) {
  .Call(.xhp__sample_partition_given_n_clusters, xhp, n_clusters)
}

#' @export
log_probability_cluster_sizes_given_n_clusters <- function(xhp, cluster_sizes) {
  .Call(.xhp__log_probability_cluster_sizes_given_n_clusters, xhp, cluster_sizes)
}

#' @export
log_probability_n_clusters <- function(xhp, n_clusters) {
  .Call(.xhp__log_probability_n_clusters, xhp, n_clusters)
}

#' @export
log_probability_partition <- function(xhp, partition) {
  .Call(.xhp__log_probability_partition, xhp, partition)
}

#' @export
log_probability_partition_given_cluster_sizes <- function(xhp, cluster_sizes) {
  .Call(.xhp__log_probability_partition_given_cluster_sizes, xhp, cluster_sizes)
}

#' @export
log_probability_partition_using_cluster_sizes <- function(xhp, cluster_sizes) {
  .Call(.xhp__log_probability_partition_using_cluster_sizes, xhp, cluster_sizes)
}

