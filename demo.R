library(gupd)

lambda <- 3
prob <- dpois(1:5, lambda)
n_items <- 11

f <- make_gupd(prob, n_items, log = FALSE)

# Two ways to use the function.
f(c(1,1,1,1,1,1,2,2,2,3,3))
f(3)

partitions <- salso::enumerate.partitions(n_items)
p <- apply(partitions, 1, f)
sum(p)

lambda <- 10
prob <- dpois(1:50, lambda)
n_items <- 1000000

k <- c(2, 5, 20, 30, 50)

system.time(f1 <- make_gupd(prob, n_items, log = TRUE))
o1 <- sapply(k, f1)
o1

sessionInfo()

