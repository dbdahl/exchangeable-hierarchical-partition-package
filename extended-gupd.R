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

y <- rsizes(n_items, n_clusters, 0.5, 1)
dsizes(y, 0.5, 1, log = TRUE)



rpartition <- function(configuration) {
  n <- sum(configuration)
  permutation <- sample(seq_len(n))
  result <- list()
  for (i in seq_along(configuration)) {
    w <- seq_len(configuration[i])
    result[[i]] <- permutation[w]
    permutation <- permutation[-w]
  }
  label <- order(sapply(result, \(x) min(x)))
  partition <- integer(n)
  for (i in seq_along(label)) {
    partition[result[[label[i]]]] <- i
  }
  partition
}

rpartition(c(14,3,2,1))

rev(sort(table(rpartition(rsizes(n_items, n_clusters, 0.1, 1)))))
rev(sort(table(rpartition(rsizes(n_items, n_clusters, 1, 0.1)))))

