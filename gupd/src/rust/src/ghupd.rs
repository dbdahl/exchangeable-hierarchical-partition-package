use roxido::*;

use ahash::AHashMap;
use rand::distr::weighted::WeightedIndex;
use rand::distr::Distribution;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use statrs::function::gamma::ln_gamma;
use std::hash::Hash;

enum ClusterSizesDistribution {
    Uniform {
        n_items: usize,
        max_n_clusters: usize,
        table: Vec<Vec<f64>>,
    },
    TiltedUniform {
        n_items: usize,
        max_n_clusters: usize,
        table: Vec<Vec<f64>>,
        tilt: f64,
    },
    TiltedBetaBinomial {
        alpha: f64,
        beta: f64,
    },
}

impl ClusterSizesDistribution {
    fn new_uniform(n_items: usize, max_n_clusters: usize) -> Self {
        let table = Self::precompute_size_configurations_table(n_items, max_n_clusters);
        Self::Uniform {
            n_items,
            max_n_clusters,
            table,
        }
    }

    fn new_beta_binomial(alpha: f64, beta: f64) -> Self {
        Self::TiltedBetaBinomial { alpha, beta }
    }

    fn update_tilt(self, tilt: f64) -> Self {
        match self {
            Self::Uniform {
                n_items,
                max_n_clusters,
                table,
            }
            | Self::TiltedUniform {
                n_items,
                max_n_clusters,
                table,
                ..
            } => {
                if tilt == 0.0 {
                    Self::Uniform {
                        n_items,
                        max_n_clusters,
                        table,
                    }
                } else {
                    Self::TiltedUniform {
                        n_items,
                        max_n_clusters,
                        table,
                        tilt,
                    }
                }
            }
            Self::TiltedBetaBinomial { .. } => {
                panic!("Not appropriate for this variant of SizeConfigurationDistribution enum.")
            }
        }
    }

    fn sample<R: Rng>(&self, n_clusters: usize, rng: &mut R) -> Result<Vec<usize>, &'static str> {
        match self {
            Self::Uniform {
                n_items,
                max_n_clusters,
                ..
            } => {
                if n_clusters > *max_n_clusters || n_clusters == 0 {
                    return Err("'n_clusters' is out of bounds");
                }
                let mut cluster_sizes = vec![0; n_clusters];
                let mut n_items = *n_items;
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
                Ok(cluster_sizes)
            }
            Self::TiltedUniform {
                n_items,
                max_n_clusters,
                tilt,
                ..
            } => {
                if n_clusters > *max_n_clusters || n_clusters == 0 {
                    return Err("'n_clusters' is out of bounds");
                }
                let mut cluster_sizes = vec![0; n_clusters];
                let mut n_items = *n_items;
                let mut min_size = 1;
                for k in (1..=n_clusters).rev() {
                    let lw = self.log_n_size_count_configurations(n_items, k, min_size);
                    let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let w = lw.iter().map(|&x| (tilt + x - max_lw).exp());
                    let weighted_index = WeightedIndex::new(w).unwrap();
                    min_size += weighted_index.sample(rng);
                    n_items -= min_size;
                    cluster_sizes[k - 1] = min_size;
                }
                Ok(cluster_sizes)
            }
            Self::TiltedBetaBinomial { .. } => Ok(vec![0]),
        }
    }

    fn log_probability(&self, cluster_sizes: &mut [usize]) -> f64 {
        match self {
            Self::Uniform {
                n_items,
                max_n_clusters,
                ..
            } => {
                let n_clusters = cluster_sizes.len();
                if n_clusters > *max_n_clusters {
                    return f64::NEG_INFINITY;
                }
                let mut n_items_working = cluster_sizes.iter().sum::<usize>();
                if n_items_working != *n_items {
                    return f64::NEG_INFINITY;
                }
                cluster_sizes.sort_unstable_by(|a, b| b.cmp(a));
                let mut min_size = 1;
                let mut sum_log_probability = 0.0;
                for k in (1..=n_clusters).rev() {
                    let lw = self.log_n_size_count_configurations(n_items_working, k, min_size);
                    let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let w = lw.iter().map(|&x| (x - max_lw).exp());
                    let Ok(weighted_index) = WeightedIndex::new(w) else {
                        return f64::NEG_INFINITY;
                    };
                    let sampled_index = cluster_sizes[k - 1] - min_size;
                    min_size += sampled_index;
                    n_items_working -= min_size;
                    sum_log_probability += weighted_index
                        .weight(sampled_index)
                        .map(|x| x.ln())
                        .unwrap_or(f64::NEG_INFINITY)
                        - weighted_index.total_weight().ln();
                }
                sum_log_probability
            }
            Self::TiltedUniform {
                n_items,
                max_n_clusters,
                tilt,
                ..
            } => {
                let n_clusters = cluster_sizes.len();
                if n_clusters > *max_n_clusters {
                    return f64::NEG_INFINITY;
                }
                let mut n_items_working = cluster_sizes.iter().sum::<usize>();
                if n_items_working != *n_items {
                    return f64::NEG_INFINITY;
                }
                cluster_sizes.sort_unstable_by(|a, b| b.cmp(a));
                let mut min_size = 1;
                let mut sum_log_probability = 0.0;
                for k in (1..=n_clusters).rev() {
                    let lw = self.log_n_size_count_configurations(n_items_working, k, min_size);
                    let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let w = lw.iter().map(|&x| (tilt + x - max_lw).exp());
                    let Ok(weighted_index) = WeightedIndex::new(w) else {
                        return f64::NEG_INFINITY;
                    };
                    let sampled_index = cluster_sizes[k - 1] - min_size;
                    min_size += sampled_index;
                    n_items_working -= min_size;
                    sum_log_probability += weighted_index
                        .weight(sampled_index)
                        .map(|x| x.ln())
                        .unwrap_or(f64::NEG_INFINITY)
                        - weighted_index.total_weight().ln();
                }
                sum_log_probability
            }
            Self::TiltedBetaBinomial { .. } => f64::NEG_INFINITY,
        }
    }

    fn precompute_size_configurations_table(
        max_extra: usize,
        max_n_clusters: usize,
    ) -> Vec<Vec<f64>> {
        let mut dp = vec![vec![f64::NEG_INFINITY; max_n_clusters + 1]; max_extra + 1];
        for j in 0..=max_n_clusters {
            dp[0][j] = 0.0;
        }
        for i in 1..=max_extra {
            for j in 1..=max_n_clusters {
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
        match self {
            Self::Uniform { table, .. } | Self::TiltedUniform { table, .. } => {
                let max_g = n / k;
                let mut results = Vec::with_capacity(max_g - w + 1);
                for g in w..=max_g {
                    let remaining = n - g;
                    if remaining < g * (k - 1) {
                        results.push(f64::NEG_INFINITY); // log(0) = -∞
                        continue;
                    }
                    let extra = remaining - g * (k - 1);
                    results.push(table[extra][k - 1]);
                }
                results
            }
            Self::TiltedBetaBinomial { .. } => {
                panic!("Not appropriate for this variant of SizeConfigurationDistribution enum.")
            }
        }
    }
}

struct GeneralizedHierarchicalUniformPartitionDistribution {
    n_items: usize,
    n_clusters_log_probability: Vec<f64>,
    cluster_sizes_distribution: ClusterSizesDistribution,
    max_n_clusters: usize,
    n_clusters_weighted_index: WeightedIndex<f64>,
    log_factorial: Vec<f64>,
}

impl GeneralizedHierarchicalUniformPartitionDistribution {
    fn new(
        n_items: usize,
        n_clusters_log_probability: &[f64],
        cluster_sizes_distribution: ClusterSizesDistribution,
    ) -> Result<Self, &'static str> {
        if n_clusters_log_probability.len() == 0 {
            return Err("There must be at least one cluster");
        }
        let max_log = n_clusters_log_probability
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let sum_exp: f64 = n_clusters_log_probability
            .iter()
            .map(|&x| (x - max_log).exp())
            .sum();
        let log_sum_exp = max_log + sum_exp.ln();
        let n_clusters_log_probability = std::iter::once(f64::NEG_INFINITY)
            .chain(n_clusters_log_probability.iter().map(|&x| x - log_sum_exp))
            .collect::<Vec<_>>();
        let max_n_clusters = n_clusters_log_probability.len() - 1;
        let n_clusters_probability = n_clusters_log_probability
            .iter()
            .map(|&x| x.exp())
            .collect::<Vec<_>>();
        // unwrap since n_clusters_probabilities are known to be okay.
        let Ok(n_clusters_weighted_index) = WeightedIndex::new(&n_clusters_probability) else {
            return Err("Invalid distribution for the number of clusters");
        };
        let mut log_factorial = Vec::with_capacity(n_items + 1);
        for i in 0..=n_items {
            log_factorial.push(ln_gamma((i as f64) + 1.0));
        }
        Ok(Self {
            n_items,
            n_clusters_log_probability,
            cluster_sizes_distribution,
            max_n_clusters,
            n_clusters_weighted_index,
            log_factorial,
        })
    }

    /// Sample a partition from the GHUP distribution.
    fn sample_partition<R: Rng>(&mut self, rng: &mut R) -> Vec<usize> {
        let n_clusters = self.sample_n_clusters(rng);
        // unwrap is okay since n_clusters came from us.
        self.sample_partition_given_n_clusters(n_clusters, rng)
            .unwrap()
    }

    /// Given a cluster size configuration, sample a partition from the GHUP distribution.
    fn sample_partition_given_n_clusters<R: Rng>(
        &mut self,
        n_clusters: usize,
        rng: &mut R,
    ) -> Result<Vec<usize>, &'static str> {
        let cluster_sizes = self.cluster_sizes_distribution.sample(n_clusters, rng)?;
        // unwrap is okay since cluster_sizes came from us.
        Ok(self
            .sample_partition_given_cluster_sizes(&cluster_sizes, rng)
            .unwrap())
    }

    /// Given a cluster size configuration, sample a partition from the GHUP distribution.
    fn sample_partition_given_cluster_sizes<R: Rng>(
        &self,
        cluster_sizes: &[usize],
        rng: &mut R,
    ) -> Result<Vec<usize>, &'static str> {
        let n_clusters = cluster_sizes.len();
        if n_clusters > self.max_n_clusters {
            return Err("'cluster_sizes' implies more clusters than specified");
        }
        let n_items = cluster_sizes.iter().sum::<usize>();
        if n_items != self.n_items {
            return Err("'cluster_sizes' implies the wrong number of items");
        }
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
        Ok(result)
    }

    /// Sample a number of clusters from the GHUP distribution.
    fn sample_n_clusters<R: Rng>(&self, rng: &mut R) -> usize {
        self.n_clusters_weighted_index.sample(rng)
    }

    /// Log probability of a partition
    fn log_probability_partition<'a, T: 'a + Eq + Hash + Copy>(&self, partition: &[T]) -> f64 {
        let mut cluster_sizes = compute_cluster_sizes(partition.iter());
        self.log_probability_partition_using_cluster_sizes(&mut cluster_sizes)
    }

    fn log_probability_partition_using_cluster_sizes(&self, cluster_sizes: &mut [usize]) -> f64 {
        let mut sum = self.log_probability_n_clusters(cluster_sizes.len());
        sum += self
            .cluster_sizes_distribution
            .log_probability(cluster_sizes);
        sum += self.log_probability_partition_given_cluster_sizes(cluster_sizes);
        sum
    }

    fn log_probability_n_clusters(&self, n_clusters: usize) -> f64 {
        if n_clusters > self.max_n_clusters {
            f64::NEG_INFINITY
        } else {
            self.n_clusters_log_probability[n_clusters]
        }
    }

    fn log_probability_partition_given_cluster_sizes(&self, cluster_sizes: &[usize]) -> f64 {
        let n_items = cluster_sizes.iter().sum::<usize>();
        if n_items != self.n_items {
            return f64::NEG_INFINITY;
        }
        let mut log_partitions = self.log_factorial[n_items];
        for &size in cluster_sizes {
            log_partitions -= self.log_factorial[size];
        }
        let mut counts = AHashMap::new();
        for &size in cluster_sizes {
            *counts.entry(size).or_insert(0) += 1;
        }
        for &count in counts.values() {
            log_partitions -= self.log_factorial[count];
        }
        -log_partitions
    }
}

fn compute_cluster_sizes<'a, I, T: 'a + Eq + Hash + Copy>(partition: I) -> Vec<usize>
where
    I: Iterator<Item = &'a T>,
{
    let mut counts = AHashMap::new();
    for &label in partition {
        *counts.entry(label).or_insert(0) += 1;
    }
    let mut cluster_sizes: Vec<_> = counts.into_values().collect();
    cluster_sizes.sort_unstable_by(|a, b| b.cmp(a));
    cluster_sizes
}

pub fn entropy_from_partition<'a, I, T: 'a + Eq + Hash + Copy>(partition: I) -> f64
where
    I: ExactSizeIterator<Item = &'a T>,
{
    let n_items = partition.len();
    entropy_from_cluster_sizes(&compute_cluster_sizes(partition), n_items)
}

pub fn entropy_from_cluster_sizes(cluster_sizes: &[usize], n_items: usize) -> f64 {
    let n_items = n_items as f64;
    cluster_sizes.iter().fold(0.0, |s, &x| {
        let p = (x as f64) / n_items;
        s - p * p.ln()
    })
}

#[roxido(module = ghupd)]
fn ghupd_new(n_items: usize, n_clusters_log_weights: &RVector, cluster_sizes_distribution: &RList) {
    let csd_name = cluster_sizes_distribution.get_by_key("name").stop();
    let csd_name = csd_name.as_scalar().stop();
    let cluster_sizes_distribution = match csd_name.str(pc) {
        "uniform" => ClusterSizesDistribution::new_uniform(n_items, n_clusters_log_weights.len()),
        "tilted_uniform" => {
            let tilt = cluster_sizes_distribution.get_by_key("tilt").stop();
            let tilt = tilt.as_scalar().stop();
            let tilt = tilt.f64();
            ClusterSizesDistribution::new_uniform(n_items, n_clusters_log_weights.len())
                .update_tilt(tilt)
        }
        "titled_beta_binomial" => {
            let get_f64 = |name: &str| -> f64 {
                let x = cluster_sizes_distribution.get_by_key(name).stop();
                let x = x.as_scalar().stop();
                x.f64()
            };
            ClusterSizesDistribution::new_beta_binomial(get_f64("alpha"), get_f64("beta"))
        }
        e => stop!("Unrecognized cluster size distribution: {}", e),
    };
    let n_clusters_log_weights = n_clusters_log_weights.to_f64(pc);
    let ghupd = GeneralizedHierarchicalUniformPartitionDistribution::new(
        n_items,
        n_clusters_log_weights.slice(),
        cluster_sizes_distribution,
    )
    .stop();
    let result = RExternalPtr::encode(ghupd, "ghupd", pc);
    result.set_class(["ghupd"].to_r(pc));
    result
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
    let partition = ghupd
        .sample_partition_given_n_clusters(n_clusters, &mut rng)
        .stop();
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
    let partition = ghupd
        .sample_partition_given_cluster_sizes(&cluster_sizes, &mut rng)
        .stop();
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
    let cluster_sizes = ghupd
        .cluster_sizes_distribution
        .sample(n_clusters, &mut rng)
        .stop();
    cluster_sizes.into_iter().map(|x| i32::try_from(x).stop())
}

#[roxido(module = ghupd)]
fn ghupd_log_probability_partition(ghupd: &mut RExternalPtr, partition: &RVector) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let partition = partition.to_i32(pc);
    let slice = partition.slice();
    ghupd.log_probability_partition(slice)
}

#[roxido(module = ghupd)]
fn ghupd_log_probability_partition_using_cluster_sizes(
    ghupd: &RExternalPtr,
    cluster_sizes: &RVector,
) {
    let ghupd = ghupd.decode_ref::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let cluster_sizes = cluster_sizes.to_i32(pc);
    let mut cluster_sizes = cluster_sizes
        .slice()
        .iter()
        .map(|&c| usize::try_from(c).stop())
        .collect::<Vec<_>>();
    ghupd.log_probability_partition_using_cluster_sizes(&mut cluster_sizes)
}

#[roxido(module = ghupd)]
fn ghupd_log_probability_n_clusters(ghupd: &RExternalPtr, n_clusters: usize) {
    let ghupd = ghupd.decode_ref::<GeneralizedHierarchicalUniformPartitionDistribution>();
    ghupd.log_probability_n_clusters(n_clusters)
}

#[roxido(module = ghupd)]
fn ghupd_log_probability_cluster_sizes_given_n_clusters(
    ghupd: &RExternalPtr,
    cluster_sizes: &RVector,
) {
    let ghupd = ghupd.decode_ref::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let cluster_sizes = cluster_sizes.to_i32(pc);
    let mut cluster_sizes = cluster_sizes
        .slice()
        .iter()
        .map(|&c| usize::try_from(c).stop())
        .collect::<Vec<_>>();
    ghupd
        .cluster_sizes_distribution
        .log_probability(&mut cluster_sizes)
}

#[roxido(module = ghupd)]
fn ghupd_log_probability_partition_given_cluster_sizes(
    ghupd: &mut RExternalPtr,
    cluster_sizes: &RVector,
) {
    let ghupd = ghupd.decode_mut::<GeneralizedHierarchicalUniformPartitionDistribution>();
    let cluster_sizes = cluster_sizes.to_i32(pc);
    let cluster_sizes = cluster_sizes
        .slice()
        .iter()
        .map(|&c| usize::try_from(c).stop())
        .collect::<Vec<_>>();
    ghupd.log_probability_partition_given_cluster_sizes(&cluster_sizes)
}
