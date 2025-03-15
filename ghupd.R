library(gupd)

n_items <- 1000

distr <- ghupd_new(10, log(c(20, 30, 10, 25, 15)), 0.0)

x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should match weights above?

distr <- ghupd_new(10, c(1, 1, 1, 1, 1), 0.0)

x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should be uniform

x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 2))), collapse="")))
x / sum(x)  # Should be uniform

x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)  # Should be uniform

x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 4))), collapse="")))
x / sum(x)  # Should be uniform


sum(sapply(seq_len(1000), \(x) length(unique(ghupd_sample_partition_given_n_clusters(distr, 4)))) != 4) # Should be zero

sum(sapply(seq_len(1000), \(x) length(unique(ghupd_sample_partition_given_n_clusters(distr, 3)))) != 3) # Should be zero

rev(sort(table(ghupd_sample_partition_given_cluster_sizes(distr, c(4, 3)))))

ghupd_sample_n_clusters(distr)
ghupd_sample_cluster_sizes_given_n_clusters(distr, 3)
table(ghupd_sample_partition_given_n_clusters(distr, 3))
length(ghupd_sample_partition_given_n_clusters(distr, 3))


rk <- function(lambda, max_n_clusters) {
  n_clusters <- Inf
  while (n_clusters > max_n_clusters) {
    n_clusters <- rpois(1, lambda)
  }
  n_clusters
}

n_clusters <- rk(10, 15)
n_clusters

rsizes(n_items, n_clusters, 0.01, 1)
rsizes(n_items, n_clusters, 0.5, 1)
rsizes(n_items, n_clusters, 1, 0.5)
rsizes(n_items, n_clusters, 1, 0.01)

x <- sapply(seq_len(1000000), \(x) {
  rsizes(n_items, n_clusters, 0.01, 1)
})
mean(apply(x, 2, \(y) identical(y, rep(125L, 8L))))

y <- rsizes(n_items, n_clusters, 0.01, 1)
dsizes(y, 0.01, 1, log = FALSE)

x <- rpartition(n_items, n_clusters, 1, 0.001)
rev(sort(table(x)))




rsizes(n_items, n_clusters, 0.01, 1)
rsizes(n_items, n_clusters, 1, 0.01)

rsizes(n_items, n_clusters, 1, 1)


for (i in 1:15) {
  print(rsizes(n_items, n_clusters, 20, 20))
}


for (i in 1:15) {
  print(rsizes(n_items, n_clusters, 50, 2))
}

