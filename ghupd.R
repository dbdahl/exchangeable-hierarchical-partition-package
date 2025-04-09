library(gupd)

# Uniform

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "uniform"))

sum(sapply(1:5, \(x) {
  exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(10 - x, x)))
})) # Should be 1


exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(9, 1)))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(8, 2)))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(7, 3)))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(6, 4)))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(5, 5)))


exp(ghupd_log_probability_n_clusters(distr, 2))

exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(9, 1)))
ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(5, 5))




distr <- ghupd_new(10, log(c(0.20, 0.30, 0.10, 0.25, 0.15)), list(method = "uniform"))

x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should match weights above

distr <- ghupd_new(10, log(c(1, 1, 1, 1, 1)), list(method = "uniform"))

x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should be uniform

# Should be an error
tryCatch(ghupd_sample_cluster_sizes_given_n_clusters(distr, 0), error = \(x) x)
ghupd_sample_cluster_sizes_given_n_clusters(distr, 5)
tryCatch(ghupd_sample_cluster_sizes_given_n_clusters(distr, 6), error = \(x) x)

x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 2))), collapse="")))
x / sum(x)  # Should be uniform

x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)  # Should be uniform

x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 4))), collapse="")))
x / sum(x)  # Should be uniform

sum(sapply(seq_len(1000), \(x) length(unique(ghupd_sample_partition_given_n_clusters(distr, 4)))) != 4) # Should be zero

sum(sapply(seq_len(1000), \(x) length(unique(ghupd_sample_partition_given_n_clusters(distr, 3)))) != 3) # Should be zero

rev(sort(table(ghupd_sample_partition_given_cluster_sizes(distr, c(6, 4)))))

# Should be an error
tryCatch(ghupd_sample_partition_given_cluster_sizes(distr, c(4, 1, 1, 1, 1, 1)), error = \(x) x)

tryCatch(ghupd_sample_partition_given_cluster_sizes(distr, c(4, 1)), error = \(x) x)

distr <- ghupd_new(10, log(c(1)), list(method = "uniform"))
tryCatch(ghupd_new(10, numeric(0), list(method = "uniform")), error = \(x) x)
tryCatch(ghupd_new(10, -Inf, list(method = "uniform")), error = \(x) x)
x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should be one

distr <- ghupd_new(10, rep(1, 5), list(method = "uniform"))
exp(ghupd_log_probability_n_clusters(distr, 0)) # Should be zero
exp(ghupd_log_probability_n_clusters(distr, 1)) # Should be 0.2
exp(ghupd_log_probability_n_clusters(distr, 5)) # Should be 0.2
exp(ghupd_log_probability_n_clusters(distr, 6)) # Should be zero

exp(-ghupd_log_probability_partition_given_cluster_sizes(distr, c(5, 3, 2))) # Should be 2520
exp(-ghupd_log_probability_partition_given_cluster_sizes(distr, c(5, 3, 1))) # Should be zero

ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(5, 3, 2)) # Should be the same
ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(6, 3, 1)) # Should be the same
ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(7, 2, 1)) # Should be the same
ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(4, 3, 3)) # Should be the same

exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(5, 5))) # Should be 0.2
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(6, 4))) # Should be 0.2
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(7, 3))) # Should be 0.2
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(8, 2))) # Should be 0.2
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(9, 1))) # Should be 0.2

distr <- ghupd_new(5, c(1, 1, 1, 1, 1), list(method = "uniform"))
cluster_sizes <- rev(sort(table(ghupd_sample_partition(distr))))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, cluster_sizes))

n_items <- 10
distr <- ghupd_new(n_items, c(1, 1, 1, 1), list(method = "uniform"))
partitions <- salso::enumerate.partitions(n_items)
x <- apply(partitions, 1, \(p) exp(ghupd_log_probability_partition(distr, p)))
sum(x == 0.0)
sum(x != 0.0)
sum(x) # Should be 1

n_items <- 10
distr <- ghupd_new(n_items, rep(1, n_items), list(method = "uniform"))
partitions <- salso::enumerate.partitions(n_items)
x <- apply(partitions, 1, \(p) exp(ghupd_log_probability_partition(distr, p)))
sum(x == 0.0)
sum(x != 0.0)
sum(x) # Should be 1


# Tilted Uniform

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = 0.0))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)  # Should be uniform

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = 0.5))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = -0.5))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = 2.0))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = -2.0))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = 20.0))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(10, log(c(1, 1, 1, 1)), list(method = "tilted_uniform", tilt = -20.0))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)


distr <- ghupd_new(10, rep(1, 5), list(method = "tilted_uniform", tilt = 0.5))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(5, 5)))
exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(9, 1)))

sum(sapply(1:5, \(x) {
  exp(ghupd_log_probability_cluster_sizes_given_n_clusters(distr, c(10 - x, x)))
})) # Should be 1



# CRP

distr <- ghupd_new(9, log(c(1, 1, 1, 1)), list(method = "crp", concentration = 2.0))
x <- table(sapply(seq_len(20000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 4))), collapse="")))
x / sum(x)


library(gourd)
n_items <- 9
concentration <- 2
crp <- CRPPartition(nItems = n_items, concentration = concentration)
x <- samplePartition(crp, 100000)
w <- apply(x, 1, \(y) length(unique(y))) == 4
mean(w)
z <- table(apply(x[w,], 1, \(y) paste0(rev(sort(table(y))), collapse="")))
z / sum(z)


# Tilted CRP

distr <- ghupd_new(9, log(c(1, 1, 1, 1)), list(method = "tilted_crp", concentration = 1.0, tilt = 0))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(9, log(c(1, 1, 1, 1)), list(method = "crp", concentration = 1.0, discount = 0.3))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(9, log(c(1, 1, 1, 1)), list(method = "tilted_crp", concentration = 1.0, tilt = 2))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)

distr <- ghupd_new(9, log(c(1, 1, 1, 1)), list(method = "tilted_crp", concentration = 1.0, tilt = 20))
x <- table(sapply(seq_len(10000), \(x) paste0(rev(sort(ghupd_sample_cluster_sizes_given_n_clusters(distr, 3))), collapse="")))
x / sum(x)




distr <- ghupd_new(1000000, rep(1, 100), list(method = "uniform"))
cluster_sizes <- rev(sort(table(ghupd_sample_partition(distr))))
ghupd_log_probability_cluster_sizes_given_n_clusters(distr, cluster_sizes)
print(system.time(sampled_partitions <- sapply(seq_len(1000), \(x) ghupd_sample_partition(distr))))
dim(sampled_partitions)
print(system.time(apply(sampled_partitions, 2, \(partition) ghupd_log_probability_partition(distr, partition))))


distr <- ghupd_new(1000000, rep(1, 100), list(method = "crp", concentration = 1.0))
cluster_sizes <- rev(sort(table(ghupd_sample_partition(distr))))
ghupd_log_probability_cluster_sizes_given_n_clusters(distr, cluster_sizes)
print(system.time(sampled_partitions <- sapply(seq_len(1000), \(x) ghupd_sample_partition(distr))))
dim(sampled_partitions)
print(system.time(apply(sampled_partitions, 2, \(partition) ghupd_log_probability_partition(distr, partition))))


engine <- function(cluster_sizes_distribution) {
  distr <- ghupd_new(50000, rep(1, 25), cluster_sizes_distribution = cluster_sizes_distribution)
  cluster_sizes <- rev(sort(table(ghupd_sample_partition(distr))))
  ghupd_log_probability_cluster_sizes_given_n_clusters(distr, cluster_sizes)
  sampled_partitions <- sapply(seq_len(1000), \(x) ghupd_sample_partition(distr))
  dim(sampled_partitions)
  apply(sampled_partitions, 2, \(partition) ghupd_log_probability_partition(distr, partition))
}

microbenchmark::microbenchmark(
  engine(list(method = "uniform")),
  engine(list(method = "tilted_uniform", tilt = 0.0000000000001)),
  engine(list(method = "crp", concentration = 1.0)),
  engine(list(method = "tilted_crp", concentration = 1.0, tilt = 0.0000000000001)),
times = 5)

