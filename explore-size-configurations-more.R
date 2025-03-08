library(gupd)

entropy(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3))
all.equal(entropy(1:10), log(10))
all.equal(entropy(rep(1, 10)), 0.0)



rbetabinom <- function(n, size, alpha, beta) {
  p <- rbeta(n, alpha, beta)
  rbinom(n, size, p)
}

maxentropy <- function(n, k) {
  result <- rep(n %/% k, k)
  remainder <- n %% k
  result[seq_len(remainder)] <- result[seq_len(remainder)] + 1
  result
}

available <- function(configuration, index) {
  previous <- if (index < length(configuration)) configuration[index + 1] else 0
  configuration[index] - previous
}

redistribute <- function(configuration, index, n_move) {
  configuration[index] <- configuration[index] - n_move
  j <- which.min(configuration[seq_len(index - 1)])
  sq <- seq(j, index - 1)
  sq <- sq[seq_len(min(n_move, length(sq)))]
  configuration[sq] <- configuration[sq] + 1
  n_move <- n_move - length(sq)
  while (n_move > 0) {
    sq <- seq_len(min(index - 1, n_move))
    configuration[sq] <- configuration[sq] + 1
    n_move <- n_move - length(sq)
  }
  configuration
}

rconfiguration <- function(n, k, alpha, beta) {
  configuration <- maxentropy(n, k) - 1
  index <- length(configuration)
  while (index > 1) {
    max_available <- available(configuration, index)
    n_move <- rbetabinom(1, max_available, alpha, beta)
    configuration <- redistribute(configuration, index, n_move)
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
