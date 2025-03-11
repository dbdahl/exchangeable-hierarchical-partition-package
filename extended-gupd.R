library(gupd)

entropy(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3))
all.equal(entropy(1:10), log(10))
all.equal(entropy(rep(1, 10)), 0.0)

# Sample from the extended gupd

n_items <- 1000

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
mean(apply(x, 2, \(y) identical(y, c(rep(91L,10L),90L))))

y <- rsizes(n_items, n_clusters, 0.01, 1)
dsizes(y, 0.01, 1, log = FALSE)

x <- rpartition(n_items, n_clusters, 1, 1)
rev(sort(table(x)))

