#!/usr/bin/env Rscript

library(gupd)

# distr <- ghupd_new(9, log(c(1, 1, 1, 1)), list(method = "crp", concentration = 1.0))
# x <- t(sapply(seq_len(10000), \(i) ghupd_sample_cluster_sizes_given_n_clusters(distr, 4)))
# w <- TRUE
# z <- table(apply(x[w,], 1, \(y) paste0(y, collapse="")))
# z / sum(z)

library(gourd)
n_items <- 10
concentration <- 2
crp <- CRPPartition(nItems = n_items, concentration = concentration)
x <- samplePartition(crp, 100000)
w <- apply(x, 1, \(y) length(unique(y))) == 3
mean(w)
z <- table(apply(x[w,], 1, \(y) paste0(rev(sort(table(y))), collapse="")))
z / sum(z)

