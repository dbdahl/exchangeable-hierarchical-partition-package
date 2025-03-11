library(gupd)

entropy(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3))
all.equal(entropy(1:10), log(10))
all.equal(entropy(rep(1, 10)), 0.0)

rgup(10000, 5, 0.001, 1)
rgup(10000, 5, 0.5, 1)
rgup(10000, 5, 0.1, 1)
rgup(10000, 5, 1, 0.5)
rgup(10000, 5, 1, 0.01)
rgup(10000, 5, 1, 0.001)


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

result <- rpartition(c(14,3,2,1))
result

size_configuration_sample(18, 5, 0.01, 1)
size_configuration_sample(18, 5, 1, 0.01)

rconfiguration(18, 5, 0.01, 1)
rconfiguration(18, 5, 1, 0.5)
rconfiguration(18, 5, 1, 0.01)

rconfiguration(180, 5, 0.01, 1)
rconfiguration(180, 5, 1, 0.5)
rconfiguration(180, 5, 1, 0.01)

rconfiguration(10000, 5, 0.001, 1)
rconfiguration(10000, 5, 0.5, 1)
rconfiguration(10000, 5, 0.1, 1)
rconfiguration(10000, 5, 1, 0.5)
rconfiguration(10000, 5, 1, 0.01)
rconfiguration(10000, 5, 1, 0.001)
