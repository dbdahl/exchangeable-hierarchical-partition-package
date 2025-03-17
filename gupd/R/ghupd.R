#' @export
ghupd_new <- function(n_items, n_clusters_log_weights, tilt) {
  .Call(.ghupd__ghupd_new, n_items, n_clusters_log_weights, tilt)
}

#' @export
print.ghupd <- function(ghupd) {
  cat("Pointer to a Generalized Hierarchical Uniform Partition Distribution (GHUPD)\n")
}

#' @export
ghupd_sample_cluster_sizes_given_n_clusters <- function(ghupd, n_clusters) {
  .Call(.ghupd__ghupd_sample_cluster_sizes_given_n_clusters, ghupd, n_clusters)
}

#' @export
ghupd_sample_n_clusters <- function(ghupd) {
  .Call(.ghupd__ghupd_sample_n_clusters, ghupd)
}

#' @export
ghupd_sample_partition <- function(ghupd) {
  .Call(.ghupd__ghupd_sample_partition, ghupd)
}

#' @export
ghupd_sample_partition_given_cluster_sizes <- function(ghupd, cluster_sizes) {
  .Call(.ghupd__ghupd_sample_partition_given_cluster_sizes, ghupd, cluster_sizes)
}

#' @export
ghupd_sample_partition_given_n_clusters <- function(ghupd, n_clusters) {
  .Call(.ghupd__ghupd_sample_partition_given_n_clusters, ghupd, n_clusters)
}

#' @export
ghupd_log_probability_cluster_sizes_given_n_clusters <- function(ghupd, cluster_sizes) {
  .Call(.ghupd__ghupd_log_probability_cluster_sizes_given_n_clusters, ghupd, cluster_sizes)
}

#' @export
ghupd_log_probability_n_clusters <- function(ghupd, n_clusters) {
  .Call(.ghupd__ghupd_log_probability_n_clusters, ghupd, n_clusters)
}

#' @export
ghupd_log_probability_partition <- function(ghupd, partition) {
  .Call(.ghupd__ghupd_log_probability_partition, ghupd, partition)
}

#' @export
ghupd_log_probability_partition_given_cluster_sizes <- function(ghupd, cluster_sizes) {
  .Call(.ghupd__ghupd_log_probability_partition_given_cluster_sizes, ghupd, cluster_sizes)
}

#' @export
ghupd_log_probability_partition_using_cluster_sizes <- function(ghupd, cluster_sizes) {
  .Call(.ghupd__ghupd_log_probability_partition_using_cluster_sizes, ghupd, cluster_sizes)
}

