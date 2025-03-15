library(gupd)

n_items <- 1000

distr <- ghupd_new(10, log(c(0.20, 0.30, 0.10, 0.25, 0.15)), 0.0)

x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should match weights above?

distr <- ghupd_new(10, c(1, 1, 1, 1, 1), 0.0)

x <- table(sapply(seq_len(10000), \(x) ghupd_sample_n_clusters(distr)))
x / sum(x)  # Should be uniform

ghupd_sample_cluster_sizes_given_n_clusters(distr, 5)
# Should be an error
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

