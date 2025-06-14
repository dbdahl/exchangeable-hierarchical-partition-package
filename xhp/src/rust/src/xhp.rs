use ahash::AHashMap;
use rand::distr::weighted::WeightedIndex;
use rand::distr::Distribution;
use rand::seq::SliceRandom;
use rand::Rng;
use roxido::*;
use statrs::function::gamma::ln_gamma;
use std::hash::Hash;

pub struct Builder {
    n_items: usize,
    max_n_clusters: usize,
    log_factorial: Vec<f64>,
}

pub struct BuilderWithNumberOfClustersDistribution {
    n_items: usize,
    max_n_clusters: usize,
    n_clusters_distribution: NumberOfClustersDistribution,
    log_factorial: Vec<f64>,
}

#[derive(Debug)]
pub struct ExchangeableHierarchicalPartitionDistribution {
    n_items: usize,
    max_n_clusters: usize,
    n_clusters_distribution: NumberOfClustersDistribution,
    cluster_sizes_distribution: ClusterSizesDistribution,
    log_factorial: Vec<f64>,
}

#[derive(Debug)]
#[allow(dead_code)]
pub enum NumberOfClustersDistribution {
    General {
        log_probability: Vec<f64>,
        weighted_index: WeightedIndex<f64>,
    },
    Crp2 {
        concentration: f64,
        discount: f64,
    },
    Binomial {
        n_trials: usize,
        probability: f64,
    },
    Poisson {
        rate: f64,
    },
    NegativeBinomial {
        n_successes: usize,
        probability: f64,
    },
}

#[derive(Debug)]
#[allow(dead_code)]
pub enum ClusterSizesDistribution {
    TiltedUniform {
        tilt: f64,
        table: Vec<Vec<f64>>,
    },
    Crp2 {
        concentration: f64,
        discount: f64,
        log_stirling: Vec<Vec<f64>>,
    },
    TiltedCrp1 {
        tilt: f64,
        log_stirling: Vec<Vec<f64>>,
    },
    TiltedBetaBinomial {
        alpha: f64,
        beta: f64,
    },
}

impl Builder {
    pub fn new(n_items: usize, max_n_clusters: usize) -> Result<Self, &'static str> {
        if n_items == 0 {
            return Err("Number of items must be at least 1.");
        };
        if max_n_clusters == 0 {
            return Err("Maximum number of clusters must be at least 1.");
        };
        if max_n_clusters > n_items {
            return Err("Maximum number of clusters cannot be greater than the number of items.");
        };
        let mut log_factorial = Vec::with_capacity(n_items + 1);
        for i in 0..=n_items {
            log_factorial.push(ln_gamma((i as f64) + 1.0));
        }
        Ok(Self {
            n_items,
            max_n_clusters,
            log_factorial,
        })
    }

    pub fn general(
        self,
        log_probability: &[f64],
    ) -> Result<BuilderWithNumberOfClustersDistribution, &'static str> {
        if self.max_n_clusters != log_probability.len() {
            return Err("Inconsistent maximum number of clusters.");
        }
        let max_log = log_probability
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let sum_exp: f64 = log_probability.iter().map(|&x| (x - max_log).exp()).sum();
        let log_sum_exp = max_log + sum_exp.ln();
        let log_probability = std::iter::once(f64::NEG_INFINITY)
            .chain(log_probability.iter().map(|&x| x - log_sum_exp))
            .collect::<Vec<_>>();
        let probability = log_probability.iter().map(|&x| x.exp()).collect::<Vec<_>>();
        let Ok(weighted_index) = WeightedIndex::new(&probability) else {
            return Err("Invalid distribution for the number of clusters.");
        };
        Ok(BuilderWithNumberOfClustersDistribution::new(
            self,
            NumberOfClustersDistribution::General {
                log_probability,
                weighted_index,
            },
        ))
    }

    #[allow(dead_code)]
    pub fn crp2(
        builder: Builder,
        concentration: f64,
        discount: f64,
    ) -> Result<BuilderWithNumberOfClustersDistribution, &'static str> {
        if discount < 0.0 || discount >= 1.0 {
            return Err("The discount parameter must be in [0,1).");
        }
        if concentration <= -discount {
            return Err("The concentration parameter must be greater than the negation of the discount parameter.");
        }
        Ok(BuilderWithNumberOfClustersDistribution::new(
            builder,
            NumberOfClustersDistribution::Crp2 {
                concentration,
                discount,
            },
        ))
    }

    #[allow(dead_code)]
    pub fn binomial(
        builder: Builder,
        n_trials: usize,
        probability: f64,
    ) -> Result<BuilderWithNumberOfClustersDistribution, &'static str> {
        if probability <= 0.0 || probability >= 1.0 {
            return Err("The probability parameter must be in [0,1].");
        }
        if n_trials <= builder.max_n_clusters {
            return Err("The number of trials must be at least the maximum number of clusters.");
        }
        Ok(BuilderWithNumberOfClustersDistribution::new(
            builder,
            NumberOfClustersDistribution::Binomial {
                n_trials,
                probability,
            },
        ))
    }

    #[allow(dead_code)]
    pub fn poisson(
        builder: Builder,
        rate: f64,
    ) -> Result<BuilderWithNumberOfClustersDistribution, &'static str> {
        if rate <= 0.0 {
            return Err("The rate parameter must be in greater than 0.0.");
        }
        Ok(BuilderWithNumberOfClustersDistribution::new(
            builder,
            NumberOfClustersDistribution::Poisson { rate },
        ))
    }

    #[allow(dead_code)]
    pub fn negative_binomial(
        builder: Builder,
        n_successes: usize,
        probability: f64,
    ) -> Result<BuilderWithNumberOfClustersDistribution, &'static str> {
        if n_successes == 0 {
            return Err("The number of success must be greater than 0.");
        }
        if probability <= 0.0 || probability >= 1.0 {
            return Err("The probability parameter must be in [0,1].");
        }
        Ok(BuilderWithNumberOfClustersDistribution::new(
            builder,
            NumberOfClustersDistribution::NegativeBinomial {
                n_successes,
                probability,
            },
        ))
    }
}

impl BuilderWithNumberOfClustersDistribution {
    pub fn new(builder: Builder, n_clusters_distribution: NumberOfClustersDistribution) -> Self {
        Self {
            n_items: builder.n_items,
            max_n_clusters: builder.max_n_clusters,
            n_clusters_distribution,
            log_factorial: builder.log_factorial,
        }
    }

    pub fn uniform(self) -> ExchangeableHierarchicalPartitionDistribution {
        let mut table = vec![vec![f64::NEG_INFINITY; self.max_n_clusters + 1]; self.n_items + 1];
        for j in 0..=self.max_n_clusters {
            table[0][j] = 0.0;
        }
        for i in 1..=self.n_items {
            for j in 1..=self.max_n_clusters {
                // Option 1: Do not use a j-sized part.
                let log_without = table[i][j - 1];
                // Option 2: Use at least one item in a part of size j.
                let log_with = if i >= j {
                    table[i - j][j]
                } else {
                    f64::NEG_INFINITY
                };
                table[i][j] = if log_without == f64::NEG_INFINITY {
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
        let dist = ClusterSizesDistribution::TiltedUniform { tilt: 0.0, table };
        ExchangeableHierarchicalPartitionDistribution {
            n_items: self.n_items,
            max_n_clusters: self.max_n_clusters,
            n_clusters_distribution: self.n_clusters_distribution,
            cluster_sizes_distribution: dist,
            log_factorial: self.log_factorial,
        }
    }

    pub fn crp1(self) -> ExchangeableHierarchicalPartitionDistribution {
        let log_stirling = Self::generate_log_stirling_table(self.n_items, self.max_n_clusters);
        let dist = ClusterSizesDistribution::TiltedCrp1 {
            tilt: 0.0,
            log_stirling,
        };
        ExchangeableHierarchicalPartitionDistribution {
            n_items: self.n_items,
            max_n_clusters: self.max_n_clusters,
            n_clusters_distribution: self.n_clusters_distribution,
            cluster_sizes_distribution: dist,
            log_factorial: self.log_factorial,
        }
    }

    pub fn crp2(
        self,
        concentration: f64,
        discount: f64,
    ) -> ExchangeableHierarchicalPartitionDistribution {
        let log_stirling = Self::generate_log_stirling_table(self.n_items, self.max_n_clusters);
        let dist = ClusterSizesDistribution::Crp2 {
            concentration,
            discount,
            log_stirling,
        };
        ExchangeableHierarchicalPartitionDistribution {
            n_items: self.n_items,
            max_n_clusters: self.max_n_clusters,
            n_clusters_distribution: self.n_clusters_distribution,
            cluster_sizes_distribution: dist,
            log_factorial: self.log_factorial,
        }
    }

    pub fn beta_binomial(
        self,
        alpha: f64,
        beta: f64,
    ) -> Result<ExchangeableHierarchicalPartitionDistribution, &'static str> {
        if alpha <= 0.0 || beta <= 0.0 {
            Err("alpha and beta must be greater than zero.")
        } else {
            let dist = ClusterSizesDistribution::TiltedBetaBinomial { alpha, beta };
            let xhp = ExchangeableHierarchicalPartitionDistribution {
                n_items: self.n_items,
                max_n_clusters: self.max_n_clusters,
                n_clusters_distribution: self.n_clusters_distribution,
                cluster_sizes_distribution: dist,
                log_factorial: self.log_factorial,
            };
            Ok(xhp)
        }
    }

    /// Generates a lookup table for log_stirling numbers with bounds on n and k.
    /// The table is a Vec<Vec<f64>> where for each n (0 <= n <= n_max)
    /// the valid k indices are 0 <= k <= min(n, k_max).
    ///
    /// The recurrence is defined as:
    /// - log_stirling(0, 0) = 0.0,
    /// - For n > 0, log_stirling(n, 0) = f64::NEG_INFINITY,
    /// - For n <= k_max, log_stirling(n, n) = 0.0,
    /// - Otherwise for 1 <= k <= min(n, k_max):
    ///   log_stirling(n, k) = log_sum_exp( log_stirling(n-1, k-1),
    ///   ln(n-1) + log_stirling(n-1, k) )
    ///
    /// # Arguments
    ///
    /// * `n_max` - maximum n value to compute.
    /// * `k_max` - maximum k value to compute (for each n, only values up to min(n, k_max) are computed).
    fn generate_log_stirling_table(n_max: usize, k_max: usize) -> Vec<Vec<f64>> {
        // Allocate table where row n has min(n, k_max)+1 entries.
        let mut table: Vec<Vec<f64>> = Vec::with_capacity(n_max + 1);
        for n in 0..=n_max {
            let num_cols = std::cmp::min(n, k_max) + 1;
            table.push(vec![f64::NEG_INFINITY; num_cols]);
        }

        // Base case: log_stirling(0, 0) = log(1) = 0.
        table[0][0] = 0.0;

        // Build the table using the recurrence.
        for n in 1..=n_max {
            let max_k = std::cmp::min(n, k_max);
            for k in 1..=max_k {
                // When n equals k (and n <= k_max), there's exactly one permutation so log(1) = 0.
                if n <= k_max && n == k {
                    table[n][k] = 0.0;
                } else {
                    // Use the recurrence:
                    // log_stirling(n, k) = log_sum_exp( log_stirling(n-1, k-1),
                    //                                   ln(n-1) + log_stirling(n-1, k) )
                    table[n][k] = Self::log_sum_exp(
                        table[n - 1][k - 1],
                        ((n - 1) as f64).ln() + table[n - 1][k],
                    );
                }
            }
        }
        table
    }

    /// Computes log(exp(a) + exp(b)) in a numerically stable way.
    fn log_sum_exp(a: f64, b: f64) -> f64 {
        if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
            return f64::NEG_INFINITY;
        }
        let m = a.max(b);
        m + ((a - m).exp() + (b - m).exp()).ln()
    }
}

impl ExchangeableHierarchicalPartitionDistribution {
    pub fn n_items(&self) -> usize {
        self.n_items
    }

    pub fn max_n_clusters(&self) -> usize {
        self.max_n_clusters
    }

    pub fn n_clusters_distribution(&self) -> &NumberOfClustersDistribution {
        &self.n_clusters_distribution
    }

    pub fn cluster_sizes_distribution(&self) -> &ClusterSizesDistribution {
        &self.cluster_sizes_distribution
    }

    pub fn cluster_sizes_distribution_mut(&mut self) -> &mut ClusterSizesDistribution {
        &mut self.cluster_sizes_distribution
    }

    /// Sample a partition from the XHP distribution.
    pub fn sample_partition<R: Rng>(&mut self, rng: &mut R) -> Vec<usize> {
        let n_clusters = self.sample_n_clusters(rng);
        // unwrap is okay since n_clusters came from us.
        self.sample_partition_given_n_clusters(n_clusters, rng)
            .unwrap()
    }

    /// Given a cluster size configuration, sample a partition from the XHP distribution.
    pub fn sample_partition_given_n_clusters<R: Rng>(
        &mut self,
        n_clusters: usize,
        rng: &mut R,
    ) -> Result<Vec<usize>, &'static str> {
        let cluster_sizes = self
            .cluster_sizes_distribution
            .sample(self, n_clusters, rng)?;
        // unwrap is okay since cluster_sizes came from us.
        Ok(self
            .sample_partition_given_cluster_sizes(&cluster_sizes, rng)
            .unwrap())
    }

    /// Given a cluster size configuration, sample a partition from the XHP distribution.
    pub fn sample_partition_given_cluster_sizes<R: Rng>(
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

    /// Sample a number of clusters from the XHP distribution.
    pub fn sample_n_clusters<R: Rng>(&self, rng: &mut R) -> usize {
        self.n_clusters_distribution.sample(rng)
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

    /// Log probability of a partition
    pub fn log_probability_partition<'a, T: 'a + Eq + Hash + Copy>(&self, partition: &[T]) -> f64 {
        let mut cluster_sizes = Self::compute_cluster_sizes(partition.iter());
        self.log_probability_partition_using_cluster_sizes(&mut cluster_sizes)
    }

    pub fn log_probability_partition_using_cluster_sizes(
        &self,
        cluster_sizes: &mut [usize],
    ) -> f64 {
        let mut sum = self.log_probability_n_clusters(cluster_sizes.len());
        sum += self
            .cluster_sizes_distribution
            .log_probability(self, cluster_sizes);
        sum += self.log_probability_partition_given_cluster_sizes(cluster_sizes);
        sum
    }

    pub fn log_probability_n_clusters(&self, n_clusters: usize) -> f64 {
        if n_clusters > self.max_n_clusters {
            f64::NEG_INFINITY
        } else {
            self.n_clusters_distribution.log_probability(n_clusters)
        }
    }

    pub fn log_probability_partition_given_cluster_sizes(&self, cluster_sizes: &[usize]) -> f64 {
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

impl NumberOfClustersDistribution {
    pub fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        match self {
            Self::General { weighted_index, .. } => weighted_index.sample(rng),
            _ => panic!("Not yet implemented."),
        }
    }

    pub fn log_probability(&self, n_clusters: usize) -> f64 {
        match self {
            Self::General {
                log_probability, ..
            } => *log_probability
                .get(n_clusters)
                .unwrap_or(&f64::NEG_INFINITY),
            _ => panic!("Not yet implemented."),
        }
    }
}

impl ClusterSizesDistribution {
    pub fn update_tilt(&mut self, x: f64) {
        match self {
            Self::TiltedUniform { tilt, .. } => *tilt = x,
            Self::TiltedCrp1 { tilt, .. } => *tilt = x,
            _ => {}
        }
    }

    pub fn sample<R: Rng>(
        &self,
        xhp: &ExchangeableHierarchicalPartitionDistribution,
        n_clusters: usize,
        rng: &mut R,
    ) -> Result<Vec<usize>, &'static str> {
        match self {
            Self::Crp2 {
                concentration,
                discount,
                ..
            } => {
                unimplemented!();
            }
            Self::TiltedUniform { tilt, .. } => {
                if n_clusters > xhp.max_n_clusters || n_clusters == 0 {
                    return Err(
                        "Number of clusters is greater than the maximum number of clusters.",
                    );
                }
                let mut cluster_sizes = vec![0; n_clusters];
                let mut n_items_working = xhp.n_items;
                let mut min_size = 1;
                for k in (1..=n_clusters).rev() {
                    let mut lw = self.log_n_size_count_configurations(n_items_working, k, min_size);
                    if *tilt != 0.0 {
                        let ideal_size_for_max_entropy = (n_items_working as f64) / (k as f64);
                        let range = min_size..(min_size + lw.len());
                        for (lw, s) in lw.iter_mut().zip(range) {
                            *lw -= *tilt * (s as f64 - ideal_size_for_max_entropy).powi(2);
                        }
                    }
                    let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let w = lw.iter().map(|&x| (x - max_lw).exp());
                    let weighted_index = WeightedIndex::new(w).unwrap();
                    min_size += weighted_index.sample(rng);
                    n_items_working -= min_size;
                    cluster_sizes[k - 1] = min_size;
                }
                Ok(cluster_sizes)
            }
            Self::TiltedCrp1 {
                tilt, log_stirling, ..
            } => {
                // Start with the full set of items and clusters.
                let mut n_items = xhp.n_items;
                let mut n_clusters = n_clusters;
                let mut result = Vec::new();
                // We'll reuse this vector for the sake of efficiency.
                let mut lw: Vec<f64> = Vec::with_capacity(n_items - n_clusters + 1);
                while n_clusters > 1 {
                    let ideal_size_for_max_entropy = (n_items as f64) / (n_clusters as f64);
                    // Each remaining cluster must get at least one item.
                    let max_s = n_items - n_clusters + 1;
                    // For each candidate first cluster size s, compute the weight.
                    for s in 1..=max_s {
                        let log_weight = xhp.log_factorial[n_items - 2]
                            - xhp.log_factorial[n_items - s]
                            + log_stirling[n_items - s][n_clusters - 1];
                        // Apply the tilt adjustment.
                        lw.push(
                            log_weight - *tilt * (s as f64 - ideal_size_for_max_entropy).powi(2),
                        );
                    }
                    // Normalize weights to avoid overflow.
                    let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let weights: Vec<f64> = lw.iter().map(|&x| (x - max_lw).exp()).collect();
                    // Sample s from the candidate sizes using the computed weights.
                    let weighted_index = WeightedIndex::new(&weights).unwrap();
                    let s_sample = 1 + weighted_index.sample(rng);
                    result.push(s_sample);
                    // Update the remaining items and cluster count.
                    n_items -= s_sample;
                    n_clusters -= 1;
                    lw.clear();
                }
                // The last cluster gets the remaining items.
                result.push(n_items);
                result.sort_unstable_by(|a, b| b.cmp(a));
                Ok(result)
            }
            Self::TiltedBetaBinomial { .. } => Err("Not yet implemented."),
        }
    }

    pub fn log_probability(
        &self,
        xhp: &ExchangeableHierarchicalPartitionDistribution,
        cluster_sizes: &mut [usize],
    ) -> f64 {
        match self {
            Self::Crp2 {
                concentration,
                discount,
                ..
            } => {
                unimplemented!();
            }
            Self::TiltedUniform { tilt, .. } => {
                let n_clusters = cluster_sizes.len();
                if n_clusters > xhp.max_n_clusters {
                    return f64::NEG_INFINITY;
                }
                let mut n_items_working = cluster_sizes.iter().sum::<usize>();
                if n_items_working != xhp.n_items {
                    return f64::NEG_INFINITY;
                }
                cluster_sizes.sort_unstable_by(|a, b| b.cmp(a));
                let mut min_size = 1;
                let mut sum_log_probability = 0.0;
                for k in (1..=n_clusters).rev() {
                    let mut lw = self.log_n_size_count_configurations(n_items_working, k, min_size);
                    if *tilt != 0.0 {
                        let ideal_size_for_max_entropy = (n_items_working as f64) / (k as f64);
                        let range = min_size..(min_size + lw.len());
                        for (lw, s) in lw.iter_mut().zip(range) {
                            *lw -= *tilt * (s as f64 - ideal_size_for_max_entropy).powi(2);
                        }
                    }
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
            Self::TiltedCrp1 {
                tilt, log_stirling, ..
            } => {
                let n_clusters = cluster_sizes.len();
                if cluster_sizes.iter().sum::<usize>() != xhp.n_items {
                    return f64::NEG_INFINITY;
                }
                // cluster_sizes.sort_unstable_by(|a, b| b.cmp(a));
                // Start with the full set of items and clusters.
                let mut n_items = xhp.n_items;
                let mut n_clusters = n_clusters;
                // We'll reuse this vector for the sake of efficiency.
                let mut lw: Vec<f64> = Vec::with_capacity(n_items - n_clusters + 1);
                let mut sum_log_probability = 0.0;
                let mut counter = 0;
                while n_clusters > 1 {
                    let ideal_size_for_max_entropy = (n_items as f64) / (n_clusters as f64);
                    // Each remaining cluster must get at least one item.
                    let max_s = n_items - n_clusters + 1;
                    // For each candidate first cluster size s, compute the weight.
                    for s in 1..=max_s {
                        let log_weight = xhp.log_factorial[n_items - 2]
                            - xhp.log_factorial[n_items - s]
                            + log_stirling[n_items - s][n_clusters - 1];
                        // Apply the tilt adjustment.
                        lw.push(
                            log_weight - *tilt * (s as f64 - ideal_size_for_max_entropy).powi(2),
                        );
                    }
                    // Normalize weights to avoid overflow.
                    let max_lw = lw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    let weights: Vec<f64> = lw.iter().map(|&x| (x - max_lw).exp()).collect();
                    // Sample s from the candidate sizes using the computed weights.
                    let Ok(weighted_index) = WeightedIndex::new(&weights) else {
                        return f64::NEG_INFINITY;
                    };
                    let s_sample = cluster_sizes[counter] - 1;
                    sum_log_probability += weighted_index
                        .weight(s_sample)
                        .map(|x| x.ln())
                        .unwrap_or(f64::NEG_INFINITY)
                        - weighted_index.total_weight().ln();
                    counter += 1;
                    n_items -= s_sample;
                    n_clusters -= 1;
                    lw.clear();

                    let s_sample = 1 + cluster_sizes[counter];

                    counter += 1;
                    // Update the remaining items and cluster count.
                    n_items -= s_sample;
                    n_clusters -= 1;
                    lw.clear();
                }
                sum_log_probability
            }
            Self::TiltedBetaBinomial { .. } => f64::NEG_INFINITY,
        }
    }

    fn log_n_size_count_configurations(&self, n: usize, k: usize, w: usize) -> Vec<f64> {
        match self {
            Self::TiltedUniform { table, .. } => {
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
            _ => {
                panic!("Not appropriate for this variant of SizeConfigurationDistribution enum.")
            }
        }
    }
}

#[roxido(module = xhp)]
fn print(xhp: &mut RExternalPtr) {
    rprintln!("Generalized Hierarchical Uniform Partition Distribution (xhp)");
    let xhp = xhp.decode_mut::<ExchangeableHierarchicalPartitionDistribution>();
    rprintln!("  + Number of items: {}", xhp.n_items());
    rprintln!("  + Maximum number of clusters: {}", xhp.max_n_clusters());
    rprint!("  + Number of clusters distribution: ");
    match &xhp.n_clusters_distribution() {
        NumberOfClustersDistribution::General { weighted_index, .. } => {
            rprintln!(
                "General({})",
                weighted_index
                    .weights()
                    .map(|x| format!("{:.3}", x))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        NumberOfClustersDistribution::Crp2 {
            concentration,
            discount,
        } => {
            rprintln!("CRP(concentration = {concentration}, discount = {discount})");
        }
        NumberOfClustersDistribution::Binomial {
            n_trials,
            probability,
        } => {
            rprintln!("Binomial(n_trials = {n_trials}, probability = {probability})");
        }
        NumberOfClustersDistribution::Poisson { rate } => {
            rprintln!("Poisson(rate = {rate})");
        }
        NumberOfClustersDistribution::NegativeBinomial {
            n_successes,
            probability,
        } => {
            rprintln!("NegativeBinomial(n_successes = {n_successes}, probability = {probability})");
        }
    }
    ()
}
