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
prob <- dpois(1:20, lambda)
n_items <- 100000

f <- make_gupd(prob, n_items, log = TRUE)
f(1)
f(2)
f(9)
f(10)
f(11)
f(20)
f(21)

