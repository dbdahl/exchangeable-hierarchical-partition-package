library(gupd)

experiment(n_items = 10, n_clusters = 3, n_samples = 12, concentration = 10.0)
experiment(n_items = 10, n_clusters = 3, n_samples = 12, concentration = 1.0)
experiment(n_items = 10, n_clusters = 3, n_samples = 12, concentration = 0.0)
experiment(n_items = 10, n_clusters = 3, n_samples = 12, concentration = -1.0)
experiment(n_items = 10, n_clusters = 3, n_samples = 12, concentration = -10.0)

