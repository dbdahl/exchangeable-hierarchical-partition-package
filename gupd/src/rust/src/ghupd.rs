use roxido::*;

use ahash::AHashMap;
use rand::distr::weighted::WeightedIndex;
use rand::distr::Distribution;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;

struct GeneralizedHierarchicalUniformPartitionDistribution {
    n_items: usize,
    n_clusters_log_probability: Vec<f64>,
    tilt: f64,
    n_clusters_weighted_index: WeightedIndex<f64>,
    size_configurations_table: Vec<Vec<f64>>,
}

impl GeneralizedHierarchicalUniformPartitionDistribution {
    fn new(n_items: usize, n_clusters_log_probability: &[f64], tilt: f64) -> Self {
        let max_log = n_clusters_log_probability
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let sum_exp: f64 = n_clusters_log_probability
            .iter()
            .map(|&x| (x - max_log).exp())
            .sum();
        let log_sum_exp = max_log + sum_exp.ln();
        let n_clusters_log_probability = n_clusters_log_probability
            .iter()
            .map(|&x| x - log_sum_exp)
            .collect::<Vec<_>>();
        let n_clusters_probability = n_clusters_log_probability
            .iter()
            .map(|&x| x.exp())
            .collect::<Vec<_>>();
        // unwrap since n_clusters_probabilities are known to be okay.
        let n_clusters_weighted_index = WeightedIndex::new(&n_clusters_probability).unwrap();
        let size_configurations_table =
            Self::precompute_size_configurations_table(n_items, n_clusters_log_probability.len());
        Self {
            n_items,
            n_clusters_log_probability,
            tilt,
            n_clusters_weighted_index,
            size_configurations_table,
        }
    }

    fn precompute_size_configurations_table(
        max_extra: usize,
        max_n_clusters: usize,
    ) -> Vec<Vec<f64>> {
        let mut dp = vec![vec![f64::NEG_INFINITY; max_n_clusters]; max_extra + 1];
        for j in 0..max_n_clusters {
            dp[0][j] = 0.0;
        }
        for i in 1..=max_extra {
            for j in 1..max_n_clusters {
                // Option 1: Do not use a j-sized part.
                let log_without = dp[i][j - 1];
                // Option 2: Use at least one item in a part of size j.
                let log_with = if i >= j {
                    dp[i - j][j]
                } else {
                    f64::NEG_INFINITY
                };
                dp[i][j] = if log_without == f64::NEG_INFINITY {
                    log_with
                } else if log_with == f64::NEG_INFINITY {
                    log_without
                } else {
                    // Log-sum-exp: log(exp(log_without) + exp(log_with))
                    let max_log = log_without.max(log_with);
                    let min_log = log_without.min(log_with);
                    max_log + (min_log - max_log).exp().ln_1p()
                };
            }
        }
        dp
    }

    fn log_n_size_count_configurations(&self, n: usize, k: usize, w: usize) -> Vec<f64> {
        let max_g = n / k;
        let mut results = Vec::with_capacity(max_g - w + 1);
        for g in w..=max_g {
            let remaining = n - g;
            if remaining < g * (k - 1) {
                results.push(f64::NEG_INFINITY); // log(0) = -∞
                continue;
            }
            let extra = remaining - g * (k - 1);
            results.push(self.size_configurations_table[extra][k - 1]);
        }
        results
    }

    /// Sample a partition from the GHUP distribution.
    fn sample_partition<R: Rng>(&mut self, rng: &mut R) -> Vec<usize> {
        let n_clusters = self.sample_n_clusters(rng);
        self.sample_partition_given_n_clusters(n_clusters, rng)
    }

    /// Given a cluster size configuration, sample a partition from the GHUP distribution.
    fn sample_partition_given_n_clusters<R: Rng>(
        &mut self,
        n_clusters: usize,
        rng: &mut R,
    ) -> Vec<usize> {
        let cluster_sizes = self.sample_cluster_sizes_given_n_clusters(n_clusters, rng);
        self.sample_partition_given_cluster_sizes(&cluster_sizes, rng)
    }

    /// Given a cluster size configuration, sample a partition from the GHUP distribution.
    fn sample_partition_given_cluster_sizes<R: Rng>(
        &self,
        cluster_sizes: &[usize],
        rng: &mut R,
    ) -> Vec<usize> {
        println!("Cluster sizes: {:?}", cluster_sizes);
        let n_clusters = cluster_sizes.len();
        let n_items = cluster_sizes.iter().sum::<usize>();
        let mut permutation: Vec<_> = (0..n_items).collect();
        permutation.shuffle(rng);
        let mut label_map = Vec::with_capacity(n_clusters);
        let mut i = 0;
        for (counter, cluster_size) in cluster_sizes.iter().enumerate() {
            // unwrap since cluster size must be at least one.
            let min = permutation[i..(i + cluster_size)].iter().min().unwrap();
            i += cluster_size;
            label_map.push((counter, *min));
        }
        label_map.sort_by(|x, y| x.1.cmp(&y.1));
        for (i, x) in label_map.iter_mut().enumerate() {
            x.1 = i;
        }
        label_map.sort_by(|x, y| x.0.cmp(&y.0));
        let mut result = vec![0; n_items];
        let mut i = 0;
        for (cluster_size, label) in cluster_sizes.iter().zip(label_map.iter().map(|x| x.1)) {
            for index in permutation[i..(i + cluster_size)].iter() {
                result[*index] = label;
            }
            i += cluster_size;
        }
        result
    }

    /// Sample a number of clusters from the GHUP distribution.
    fn sample_n_clusters<R: Rng>(&self, rng: &mut R) -> usize {
        self.n_clusters_weighted_index.sample(rng) + 1
    }

    /// Given a number of clusters, sample a cluster size configuration from the GHUP distribution.
    fn sample_cluster_sizes_given_n_clusters<R: Rng>(
        &mut self,
        n_clusters: usize,
        rng: &mut R,
    ) -> Vec<usize> {
        let mut cluster_sizes = vec![0; n_clusters];
        let mut n_items = self.n_items;
        let mut min_size = 1;
        for k in (1..=n_clusters).rev() {
            let lw = self.log_n_size_count_configurations(n_items, k, min_size);
            let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let w = lw.iter().map(|&x| (x - max_lw).exp());
            let weighted_index = WeightedIndex::new(w).unwrap();
            min_size += weighted_index.sample(rng);
            n_items -= min_size;
            cluster_sizes[k - 1] = min_size;
        }
        cluster_sizes
    }

    /// Log probability of a partition
    fn log_probabilty_partition(&self, partition: &[usize]) -> f64 {
        let cluster_sizes = compute_cluster_sizes(partition);
        self.log_probabilty_partition_using_cluster_sizes(&cluster_sizes)
    }

    fn log_probabilty_partition_using_cluster_sizes(&self, cluster_sizes: &[usize]) -> f64 {
        let mut sum = self.log_probability_n_clusters(cluster_sizes.len());
        sum += self.log_probability_cluster_sizes_given_n_clusters(&cluster_sizes);
        sum += self.log_probability_partition_given_cluster_sizes(&cluster_sizes);
        sum
    }

    fn log_probability_n_clusters(&self, n_clusters: usize) -> f64 {
        self.n_clusters_log_probability[n_clusters]
    }

    fn log_probability_cluster_sizes_given_n_clusters(&self, cluster_sizes: &[usize]) -> f64 {
        todo!()
    }

    fn log_probability_partition_given_cluster_sizes(&self, cluster_sizes: &[usize]) -> f64 {
        todo!()
    }
}

pub struct SizesCounter {
    cache: AHashMap<(usize, usize), f64>,
}

impl SizesCounter {
    pub fn new() -> Self {
        Self {
            cache: AHashMap::new(),
        }
    }

    pub fn log_n_size_configurations_vector(
        &mut self,
        n_items: usize,
        n_clusters: usize,
        min_size: usize,
    ) -> Vec<f64> {
        let max_g = n_items / n_clusters;
        let mut results = Vec::with_capacity(max_g - min_size + 1);

        for g in min_size..=max_g {
            let remaining = n_items - g;

            if remaining < g * (n_clusters - 1) {
                results.push(f64::NEG_INFINITY); // log(0) = -∞
                continue;
            }

            let extra = remaining - g * (n_clusters - 1);
            let log_count = self.log_n_size_configurations(extra, n_clusters - 1);
            results.push(log_count);
        }
        results
    }

    /// Calculate or retrieve from cache the log of number of partitions
    pub fn log_n_size_configurations(&mut self, n_items: usize, n_clusters: usize) -> f64 {
        let key = (n_items, n_clusters);

        // Return cached result if available
        if let Some(&result) = self.cache.get(&key) {
            return result;
        }

        // Otherwise calculate and cache the result
        let result = self.calculate_log_partitions(n_items, n_clusters);
        self.cache.insert(key, result);

        result
    }

    /// Core calculation function with memory optimization
    fn calculate_log_partitions(&self, n: usize, k: usize) -> f64 {
        if n == 0 {
            return 0.0; // log(1) = 0
        }
        if k == 0 {
            return f64::NEG_INFINITY; // log(0) = -∞
        }

        // Two-row approach for memory efficiency
        let mut prev_row = vec![0.0];
        prev_row.resize(k + 1, 0.0);

        let mut curr_row = vec![f64::NEG_INFINITY];
        curr_row.resize(k + 1, f64::NEG_INFINITY);

        for i in 1..=n {
            for j in 1..=k {
                // Case 1: Don't use part j
                let log_without = curr_row[j - 1];

                // Case 2: Use at least one item in part j
                let log_with = if i >= j {
                    prev_row[j]
                } else {
                    f64::NEG_INFINITY // Not possible
                };

                // Combine using log-sum-exp trick
                curr_row[j] = self.log_sum_exp(log_without, log_with);
            }

            // Swap rows for next iteration
            std::mem::swap(&mut prev_row, &mut curr_row);
            curr_row[0] = f64::NEG_INFINITY;
            for j in 1..=k {
                curr_row[j] = f64::NEG_INFINITY;
            }
        }

        prev_row[k]
    }

    /// Helper function to add logarithms using the log-sum-exp trick
    fn log_sum_exp(&self, log_a: f64, log_b: f64) -> f64 {
        if log_a == f64::NEG_INFINITY {
            return log_b;
        }
        if log_b == f64::NEG_INFINITY {
            return log_a;
        }

        let max_log = log_a.max(log_b);
        let min_log = log_a.min(log_b);
        max_log + (min_log - max_log).exp().ln_1p()
    }

    /// Get current cache size
    pub fn cache_size(&self) -> usize {
        self.cache.len()
    }

    /// Clear the cache to free memory
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }
}

fn compute_cluster_sizes(partition: &[usize]) -> Vec<usize> {
    let mut counts = AHashMap::new();
    for &label in partition {
        *counts.entry(label).or_insert(0) += 1;
    }
    let mut cluster_sizes: Vec<_> = counts.into_values().collect();
    cluster_sizes.sort_unstable_by(|a, b| b.cmp(a));
    cluster_sizes
}

#[roxido(module = ghupd)]
fn ghupd_new(n_items: usize, n_clusters_log_weights: &RVector, tilt: f64) {
    let n_clusters_log_weights = n_clusters_log_weights.to_f64(pc);
    let ghupd = GeneralizedHierarchicalUniformPartitionDistribution::new(
        n_items,
        n_clusters_log_weights.slice(),
        tilt,
    );
    RExternalPtr::encode(ghupd, "ghupd", pc)
}

#[roxido(module = ghupd)]
fn ghupd_sample_partition(ghupd: &mut RExternalPtr) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let partition = ghupd.sample_partition(&mut rng);
    partition.into_iter().map(|x| i32::try_from(x).stop() + 1)
}

#[roxido(module = ghupd)]
fn ghupd_sample_partition_given_n_clusters(ghupd: &mut RExternalPtr, n_clusters: usize) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let partition = ghupd.sample_partition_given_n_clusters(n_clusters, &mut rng);
    partition.into_iter().map(|x| i32::try_from(x).stop() + 1)
}

#[roxido(module = ghupd)]
fn ghupd_sample_partition_given_cluster_sizes(ghupd: &mut RExternalPtr, cluster_sizes: &RVector) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let cluster_sizes = cluster_sizes.to_i32(pc);
    let cluster_sizes = cluster_sizes
        .slice()
        .iter()
        .map(|&x| usize::try_from(x).stop())
        .collect::<Vec<_>>();
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let partition = ghupd.sample_partition_given_cluster_sizes(&cluster_sizes, &mut rng);
    partition.into_iter().map(|x| i32::try_from(x).stop() + 1)
}

#[roxido(module = ghupd)]
fn ghupd_sample_n_clusters(ghupd: &mut RExternalPtr) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let cluster_sizes = ghupd.sample_n_clusters(&mut rng);
    i32::try_from(cluster_sizes).stop()
}

#[roxido(module = ghupd)]
fn ghupd_sample_cluster_sizes_given_n_clusters(ghupd: &mut RExternalPtr, n_clusters: usize) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let cluster_sizes = ghupd.sample_cluster_sizes_given_n_clusters(n_clusters, &mut rng);
    cluster_sizes.into_iter().map(|x| i32::try_from(x).stop())
}
