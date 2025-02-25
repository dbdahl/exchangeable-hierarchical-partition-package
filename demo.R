library(gupd)

prob <- c(2, 100, 2, 1, 1)
n_items <- 11

f <- make_gupd(prob, n_items)

# Two ways to use the function.
f(c(1,1,1,1,1,1,2,2,2,3,3))
f(3)

partitions <- salso::enumerate.partitions(n_items)
p <- apply(partitions, 1, f)
sum(p)


prob <- dpois(1:20, lambda)
n_items <- 82

f <- make_gupd(prob, n_items)
f(20)
f(21)

