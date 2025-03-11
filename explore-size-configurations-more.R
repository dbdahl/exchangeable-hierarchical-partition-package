library(gupd)

entropy(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3))
all.equal(entropy(1:10), log(10))
all.equal(entropy(rep(1, 10)), 0.0)

rsizes(10000, 5, 0.01, 1)
rsizes(10000, 5, 0.5, 1)
rsizes(10000, 5, 1, 0.5)
rsizes(10000, 5, 1, 0.01)

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

rev(sort(table(rpartition(rgup(100, 5, 0.1, 1)))))
rev(sort(table(rpartition(rgup(100, 5, 1, 0.1)))))


