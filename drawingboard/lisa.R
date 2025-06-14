gini <- function(s) {
  k <- length(s)
  p <- s / sum(s)
  1 - sum(p^2)
}
f <- gini
f <- xhp::entropy_from_cluster_sizes

rescale <- function(x) {
  x_min <- min(x)
  x_max <- max(x)
  r <- (x - x_min) / (x_max - x_min)
  if (length(r) == 1) 1.0 else r
}


n_items <- 10
k <- 3
target <- 0.57
temperature <- 10

state <- integer(k)
index <- k
min_n_items <- 1

n_items_remaining <- n_items
k_remaining <- k
while (k_remaining > 1) {
  candidates <- lapply(seq(min_n_items, floor(n_items_remaining/k_remaining)), \(x) c(n_items_remaining - x, x))
  y <- sapply(candidates, \(x) f(x))
  z <- rescale(y)
  d <- -abs(z - target)
  e <- exp(temperature * d)
  e <- e/sum(e)
  cat("---\n")
  print(k_remaining)
  print(candidates)
  print(e)
  w <- sample(seq_along(candidates), 1, prob = e)
  count <- candidates[[w]][2]
  state[k_remaining] <- count
  min_n_items <- max(min_n_items, count)
  n_items_remaining <- n_items_remaining - count
  k_remaining <- k_remaining - 1
}
state[1] <- n_items_remaining
state
f(state)
sum(state)

source("allpossible.R")
all <- enumerate_partition_configs(n_items, k) 
all <- cbind(all, apply(all, 1, \(x) f(x)))
all <- cbind(all, rescale(all[, ncol(all)]))
all <- cbind(all, apply(all[, seq_len(k)], 1, \(x) identical(x, state)))
all <- cbind(all, 0)
all[which.min(abs(all[, ncol(all) - 2] - target)), ncol(all)] <- 1
all


