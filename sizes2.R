rbetabinomial <- function(n, alpha, beta) {
  p <- rbeta(1, alpha, beta)
  p <- 1/4
  rbinom(1, n, p)
}

rsizes <- function(n_items, n_clusters, alpha, beta, method = 2) {
  sizes <- integer(n_clusters)
  n <- n_items
  min <- 1
  if (method == 1) {
    p <- rbeta(1, alpha, beta)
    cat("p: ", p, "\n", sep="")
  }
  for (i in rev(seq_len(n_clusters)[-1])) {
    max <- n %/% i
    sizes[i] <- min + if (method == 1) rbinom(1, max - min, p) else rbetabinomial(max - min, alpha, beta)
    min <- sizes[i]
    n <- n - min
  }
  sizes[1] <- n
  sizes
}

# debug(rsizes)
rsizes(100, 8, 0.01, 1)
rsizes(100, 8, 1, 0.01)
rsizes(100, 8, 0.001, 0.01)
rsizes(100, 8, 0.05, 0.1)
rsizes(100, 8, 1, 1)


rsizes(100, 8, 1, 0.00001)
rsizes(100, 8, 1, 1000000)


rsizes(1000, 8, 10, 1)
rsizes(1000, 8, 10, 10)


x <- sapply(seq_len(100000), \(x) { rsizes(8, 3, 1e100, 1e100) })
table(apply(x, 2, \(y) paste0(y, collapse="")))




