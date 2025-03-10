library(gupd)

entropy(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3))
all.equal(entropy(1:10), log(10))
all.equal(entropy(rep(1, 10)), 0.0)

x <- new_max_entropy(10, 3)
x <- new_from_vec(c(15,9,4,2,1))
x <- new_from_vec(c(15,9,4,2,0))
x <- new_from_vec(c(15,9,4,1,2))

x <- new_from_vec(c(15,9,4,30,2))

size_configuration_to_r(x)
size_configuration_available(x, 6)
size_configuration_available(x, 5)
size_configuration_available(x, 4)
size_configuration_available(x, 3)
size_configuration_available(x, 2)
size_configuration_available(x, 1)

rbetabinom <- function(n, size, alpha, beta) {
  p <- rbeta(n, alpha, beta)
  rbinom(n, size, p)
}

rconfiguration <- function(n, k, alpha, beta) {
  configuration <- size_configuration_to_r(new_max_entropy(n, k))
  index <- length(configuration)
  while (index > 1) {
    max_available <- size_configuration_available(new_from_vec(configuration+1), index)
    n_move <- rbetabinom(1, max_available, alpha, beta)
    # configuration <- redistribute(configuration, index, n_move)
    ptr <- new_from_vec(configuration + 1)
    size_configuration_redistribute(ptr, index, n_move)
    configuration <- size_configuration_to_r(ptr)
    index <- index - 1
  }
  configuration <- configuration + rep(1, length(configuration))
  configuration
}

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

rconfiguration(18, 5, 0.01, 1)
rconfiguration(18, 5, 1, 0.5)
rconfiguration(18, 5, 1, 0.01)

rconfiguration(180, 5, 0.01, 1)
rconfiguration(180, 5, 1, 0.5)
rconfiguration(180, 5, 1, 0.01)

rconfiguration(10000, 5, 0.01, 1)
rconfiguration(10000, 5, 0.001, 1)
rconfiguration(10000, 5, 0.01, 1)
rconfiguration(10000, 5, 1, 0.5)
rconfiguration(10000, 5, 1, 0.01)
rconfiguration(10000, 5, 1, 0.001)
