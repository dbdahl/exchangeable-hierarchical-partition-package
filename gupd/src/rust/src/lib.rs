// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

mod ghupd;

use ahash::AHashMap;
use rand::{Rng, SeedableRng};
use rand_distr::{Beta, BetaError, Binomial};
use rand_pcg::Pcg64Mcg;
use roxido::*;
use statrs::function::gamma::ln_gamma;
use std::f64;

#[roxido]
fn make_gupd(prob: &RVector, n_items: i32) {
    let n_items = usize::try_from(n_items).stop();
    if n_items < prob.len() {
        stop!("'n_items' should be at least the length of 'prob'.")
    }
    let prob = prob.to_f64(pc);
    let mut prob = prob.slice().to_vec();
    if prob.iter().any(|&x| x < 0.0) {
        stop!("'prob' must have all positive values")
    }
    if prob.iter().any(|&x| x.is_infinite()) {
        stop!("'prob' must have finite elements")
    }
    let sum: f64 = prob.iter().sum();
    for p in prob.iter_mut() {
        *p /= sum;
    }
    let k_max = prob.len();
    let rval = RVector::from_value(0.0, k_max, pc);
    let slice = rval.slice_mut();
    // Work entirely in the log domain.
    // The recurrence:
    //   S(n,k) = k * S(n-1,k) + S(n-1,k-1)
    // becomes:
    //   L(n,k) = log_sum_exp( ln(k) + L(n-1,k), L(n-1,k-1) )
    // Use two rows (instead of one) to allow the potential optimizations.
    let mut prev = vec![-f64::INFINITY; k_max + 1];
    let mut curr = vec![-f64::INFINITY; k_max + 1];
    // Base case: S(0,0)=1 so L(0,0)=0.
    // For k > 0, S(0,k)=0 so we set L(0,k) = negative infinity.
    prev[0] = 0.0;
    for n in 1..=n_items {
        let limit = k_max.min(n);
        curr[0] = -f64::INFINITY;
        for k in 1..=limit {
            curr[k] = log_sum_exp((k as f64).ln() + prev[k], prev[k - 1]);
        }
        // For indices k > limit, curr[k] remains negative infinity.
        // Swap the rows so that prev becomes the row for n.
        std::mem::swap(&mut prev, &mut curr);
    }
    // Note that prev holds L(n_items, k) for k = 0..k_max.
    // Skip index 0 and copy k = 1..k_max.
    for (r, (v, p)) in slice.iter_mut().zip(prev.iter().skip(1).zip(prob.iter())) {
        *r = p.ln() - v;
    }
    rval
}

/// Returns ln(exp(a) + exp(b)) in a numerically stable way.
fn log_sum_exp(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if m.is_infinite() && m.is_sign_negative() {
        -f64::INFINITY // Both a and b are negative infinity.
    } else {
        m + ((a - m).exp() + (b - m).exp()).ln()
    }
}

use rand::distr::weighted::WeightedIndex;
use rand::prelude::*;

#[roxido]
fn experiment(n_items: usize, n_clusters: usize, n_samples: usize, concentration: f64) {
    let result = RMatrix::from_value(0, n_clusters, n_samples, pc);
    let slice = result.slice_mut();
    let mut rng = rand::rng();
    let mut counts = vec![0_i32; n_clusters];
    for i in 0..n_samples {
        counts.fill(1);
        let mut r = n_items - n_clusters;
        while r > 0 {
            let weights: Vec<_> = counts
                .iter()
                .map(|&x| (x as f64).powf(concentration))
                .collect();
            let dist = WeightedIndex::new(&weights).stop();
            let index = dist.sample(&mut rng);
            counts[index] += 1;
            r -= 1;
        }
        counts.sort();
        slice[i * n_clusters..(i + 1) * n_clusters].copy_from_slice(&counts);
    }
    result
}

#[roxido]
fn entropy(partition: &RVector) {
    let partition = partition.to_i32(pc);
    let slice = partition.slice();
    let n_items = slice.len() as f64;
    let mut counts = AHashMap::new();
    for &num in slice {
        *counts.entry(num).or_insert(0) += 1;
    }
    counts.values().fold(0.0, |s, &x| {
        let p = (x as f64) / n_items;
        s - p * p.ln()
    })
}

#[roxido]
fn rsizes(n_items: i32, n_clusters: i32, alpha: f64, beta: f64) {
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let sc = SizeConfiguration::sample_size_configuration(
        u64::try_from(n_items).stop(),
        u64::try_from(n_clusters).stop(),
        alpha,
        beta,
        &mut rng,
    )
    .stop();
    let result = RVector::<i32>::new(usize::try_from(n_clusters).stop(), pc);
    let slice = result.slice_mut();
    for (x, y) in slice.iter_mut().zip(sc.x.iter()) {
        *x = i32::try_from(*y + 1).stop()
    }
    result
}

#[roxido]
fn dsizes(x: &RVector, alpha: f64, beta: f64, log: bool) {
    let mut y = Vec::with_capacity(x.len());
    if let Ok(x) = x.as_f64() {
        for z in x.slice() {
            y.push(*z as u64);
        }
    } else if let Ok(x) = x.as_i32() {
        for z in x.slice() {
            y.push(u64::try_from(*z).stop());
        }
    } else {
        stop!("Unsupported vector type");
    }
    let sc = SizeConfiguration::new_from_vec(y).stop();
    let lp = sc.log_probability(alpha, beta);
    if log {
        lp
    } else {
        lp.exp()
    }
}

#[roxido]
fn rpartition(n_items: i32, n_clusters: i32, alpha: f64, beta: f64) {
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let sc = SizeConfiguration::sample_size_configuration(
        u64::try_from(n_items).stop(),
        u64::try_from(n_clusters).stop(),
        alpha,
        beta,
        &mut rng,
    )
    .stop();
    let x = sc.sample_partition(&mut rng);
    let result = RVector::<i32>::new(usize::try_from(n_items).stop(), pc);
    let slice = result.slice_mut();
    for (z, y) in slice.iter_mut().zip(x.iter()) {
        *z = i32::try_from(*y + 1).stop()
    }
    result
}

#[derive(Debug)]
pub struct SizeConfiguration {
    x: Vec<u64>,
}

fn sample_beta_binomial<R: Rng + ?Sized>(
    n_items: u64,
    alpha: f64,
    beta: f64,
    rng: &mut R,
) -> Result<u64, BetaError> {
    let beta_dist = Beta::new(alpha, beta)?;
    let p: f64 = beta_dist.sample(rng);
    let binomial_dist =
        Binomial::new(n_items, p).expect("Beta distribution provided an invalid probability");
    Ok(binomial_dist.sample(rng))
}

pub fn log_pmf_beta_binomial(x: u64, n_items: u64, alpha: f64, beta: f64) -> f64 {
    if x > n_items {
        return f64::NEG_INFINITY;
    }
    let n = n_items as f64;
    let x_f = x as f64;
    // Compute log of the binomial coefficient: ln(n choose x)
    let log_comb = ln_gamma(n + 1.0) - ln_gamma(x_f + 1.0) - ln_gamma(n - x_f + 1.0);
    // Compute log beta function for numerator: ln(B(alpha + x, beta + n - x))
    let log_beta_num =
        ln_gamma(alpha + x_f) + ln_gamma(beta + (n - x_f)) - ln_gamma(alpha + beta + n);
    // Compute log beta function for denominator: ln(B(alpha, beta))
    let log_beta_den = ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(alpha + beta);
    log_comb + log_beta_num - log_beta_den
}

impl SizeConfiguration {
    pub fn new_from_vec(mut x: Vec<u64>) -> Result<Self, &'static str> {
        x.sort_unstable_by(|a, b| b.cmp(a));
        if x[x.len() - 1] == 0 {
            Err("All clusters must have at least one element")
        } else {
            for y in x.iter_mut() {
                *y -= 1;
            }
            Ok(Self { x })
        }
    }

    pub fn new_max_entropy(n_items: u64, n_clusters: u64) -> Result<Self, &'static str> {
        if n_clusters > n_items {
            return Err("'n_items' must be greater than 'n_clusters'");
        }
        if n_clusters == 0 {
            return Err("'n_clusters' must be greater than 0");
        }
        let min_size = n_items / n_clusters - 1;
        let mut x = vec![min_size; usize::try_from(n_clusters).expect("u64 doesn't fit in usize")];
        let remainder = n_items % n_clusters;
        for y in x[0..usize::try_from(remainder).expect("u64 doesn't fit in usize")].iter_mut() {
            *y += 1;
        }
        Ok(Self { x })
    }

    pub fn sample_size_configuration<R: Rng + ?Sized>(
        n_items: u64,
        n_clusters: u64,
        alpha: f64,
        beta: f64,
        rng: &mut R,
    ) -> Result<Self, &'static str> {
        let mut sc = Self::new_max_entropy(n_items, n_clusters)?;
        let mut index = sc.x.len() - 1;
        while index > 0 {
            let available = sc.available(index);
            let n = sample_beta_binomial(available, alpha, beta, rng)
                .map_err(|_| "Invalid beta parameter")?;
            sc.redistribute(index, n);
            index -= 1;
        }
        Ok(sc)
    }

    pub fn sample_partition<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec<u64> {
        let n_clusters = self.x.len();
        let n_items = usize::try_from(self.x.iter().sum::<u64>())
            .expect("u64 doesn't fit in usize")
            + n_clusters;
        let mut permutation: Vec<_> = (0..n_items).collect();
        permutation.shuffle(rng);
        let mut label_map = Vec::with_capacity(n_clusters);
        let mut i = 0;
        for (counter, cluster_size) in self.x.iter().enumerate() {
            let cluster_size =
                usize::try_from(*cluster_size).expect("u64 doesn't fit in usize") + 1;
            let min = permutation[i..(i + cluster_size)].iter().min().unwrap(); // Cluster size must be at least one.
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
        for (cluster_size, label) in self.x.iter().zip(label_map.iter().map(|x| x.1)) {
            let cluster_size =
                usize::try_from(*cluster_size).expect("u64 doesn't fit in usize") + 1;
            let label = u64::try_from(label).expect("usize doesn't fit in u64");
            for index in permutation[i..(i + cluster_size)].iter() {
                result[*index] = label;
            }
            i += cluster_size;
        }
        result
    }

    fn log_probability(&self, alpha: f64, beta: f64) -> f64 {
        let n_clusters = u64::try_from(self.x.len()).expect("usize doesn't fit in u64");
        let n_items = u64::try_from(self.x.iter().sum::<u64>()).expect("usize doesn't fit in u64")
            + n_clusters;
        let mut sc = Self::new_max_entropy(n_items, n_clusters).unwrap();
        let mut result = 0.0;
        let mut index = sc.x.len() - 1;
        while index > 0 {
            let available = sc.available(index);
            let n = sc.x[index] - self.x[index];
            result += log_pmf_beta_binomial(n, available, alpha, beta);
            sc.redistribute(index, n);
            index -= 1;
        }
        assert!(sc.x == self.x);
        result
    }

    fn log_partition_count_unordered(&self) -> f64 {
        let n_items: u64 = self.x.iter().sum();
        let log_n_fact = ln_gamma((n_items as f64) + 1.0);
        let log_denom: f64 = self
            .x
            .iter()
            .map(|&size| ln_gamma((size as f64) + 1.0))
            .sum();
        let mut freq: AHashMap<u64, u64> = AHashMap::new();
        for &size in self.x.iter() {
            *freq.entry(size).or_insert(0) += 1;
        }
        let log_mult: f64 = freq
            .values()
            .map(|&count| ln_gamma((count as f64) + 1.0))
            .sum();
        log_n_fact - log_denom - log_mult
    }

    fn available(&self, index: usize) -> u64 {
        match index {
            0 => 0,
            i if i >= self.x.len() => 0,
            i if i == self.x.len() - 1 => self.x[i],
            i => self.x[i] - self.x[i + 1],
        }
    }

    fn redistribute(&mut self, index: usize, mut n: u64) -> bool {
        if n == 0 {
            return true;
        }
        if index == 0 {
            return false;
        }
        self.x[index] -= n;
        let v0 = self.x[0];
        let mut i = index - 1;
        while self.x[i] != v0 {
            i -= 1;
        }
        i += 1;
        while i < index && n > 0 {
            self.x[i] += 1;
            i += 1;
            n -= 1;
        }
        if n == 0 {
            return true;
        }
        let whole = n / u64::try_from(index).expect("usize doesn't fit in u64");
        if whole > 0 {
            for j in 0..index {
                self.x[j] += whole;
            }
            n -= u64::try_from(index).expect("usize doesn't fit in u64") * whole;
        }
        for j in 0..usize::try_from(n).expect("u64 doesn't fit in usize") {
            self.x[j] += 1;
        }
        true
    }
}

fn ln_add_exp(a: f64, b: f64) -> f64 {
    // Computes ln(exp(a) + exp(b)) in a numerically stable way.
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let max_val = a.max(b);
    max_val + ((a - max_val).exp() + (b - max_val).exp()).ln()
}

fn ln_sub_exp(a: f64, b: f64) -> f64 {
    // Computes ln(exp(a) - exp(b)) assuming a >= b,
    // with special handling when b is -∞ (i.e. exp(b) = 0).
    if b == f64::NEG_INFINITY {
        return a;
    }
    let max_val = a.max(b);
    max_val + ((a - max_val).exp() - (b - max_val).exp()).ln()
}

/// Compute the natural logarithm of the number of valid partitions
/// for `n` items into `k` clusters, considering only configurations
/// where the last cluster size g is in the range w..=n/k.
/// Returns a vector of (g, ln(count)) pairs.
fn count_partitions_ln(n: usize, k: usize, w: usize) -> Vec<(usize, f64)> {
    // We'll work in log-space where ln(0) is represented by -∞.
    // We need to store dp values for each total number of items m (0..=n)
    // for the current number of clusters.
    let mut dp_prev = vec![f64::NEG_INFINITY; n + 1];
    let mut dp_curr = vec![f64::NEG_INFINITY; n + 1];
    // We'll reuse this array for the cumulative (prefix) sums.
    let mut prefix_sum_ln = vec![f64::NEG_INFINITY; n + 1];

    // Base case: With one cluster, there's exactly one way (ln(1)=0)
    for m in 0..=n {
        dp_prev[m] = 0.0;
    }

    // Build up the solution for clusters = 2 to k.
    for clusters in 2..=k {
        // Compute prefix sums for dp_prev in log-space.
        prefix_sum_ln[0] = dp_prev[0];
        for m in 1..=n {
            prefix_sum_ln[m] = ln_add_exp(prefix_sum_ln[m - 1], dp_prev[m]);
        }
        // Fill dp_curr for current number of clusters.
        // For each m (total items), dp_curr[m] represents ln(# partitions)
        // into 'clusters' clusters.
        for m in 0..=n {
            if m < clusters {
                // Not enough items to have 'clusters' clusters.
                dp_curr[m] = f64::NEG_INFINITY;
            } else {
                // Let min_size be the minimum possible size for a cluster,
                // which is given by integer division (ensuring non-increasing order).
                let min_size = m / clusters;
                // We use the prefix sums to compute the range sum in log-space.
                dp_curr[m] = ln_sub_exp(
                    prefix_sum_ln[m - clusters + 1],
                    if min_size > 0 {
                        prefix_sum_ln[min_size - 1]
                    } else {
                        f64::NEG_INFINITY
                    },
                );
            }
        }
        // Move the current dp row into dp_prev for the next iteration.
        dp_prev.copy_from_slice(&dp_curr);
    }

    // Now dp_prev is dp[k]. We want to return ln(count) for configurations
    // where the last cluster has size g, for g in w..=n/k.
    let max_g = n / k;
    let mut results = Vec::new();
    for g in w..=max_g {
        // dp_prev[n - g] holds the log count for partitions of n items
        // when the last cluster (of size g) is removed.
        results.push((g, dp_prev[n - g]));
    }
    results
}

/// Corrected function to count the logarithm of partitions with non-increasing order
fn count_configurations_log_iterative(n: usize, k: usize, w: usize) -> Vec<f64> {
    let max_g = n / k;
    let mut results = Vec::with_capacity(max_g - w + 1);

    for g in w..=max_g {
        // When g is the size of the last cluster
        let remaining = n - g;

        // For valid configurations, we need to distribute remaining items
        // into k-1 clusters, where each cluster has at least g items
        if remaining < g * (k - 1) {
            results.push(f64::NEG_INFINITY); // log(0) = -∞
            continue;
        }

        // Now we need to find partitions of (remaining) with k-1 parts
        // where each part is at least g, and parts are in non-increasing order

        // This is equivalent to finding partitions of (remaining - g*(k-1))
        // into at most k-1 parts, with no minimum size constraint
        let extra = remaining - g * (k - 1);

        // Calculate log of number of partitions of 'extra' into at most (k-1) parts
        let log_count = log_partitions_at_most_parts(extra, k - 1);
        results.push(log_count);
    }

    results
}

/// Calculate logarithm of number of partitions of n into at most k parts
/// using dynamic programming
fn log_partitions_at_most_parts(n: usize, k: usize) -> f64 {
    if n == 0 {
        return 0.0; // log(1) = 0
    }
    if k == 0 {
        return f64::NEG_INFINITY; // log(0) = -∞
    }

    // dp[i][j] = log of number of partitions of i into at most j parts
    let mut dp = vec![vec![f64::NEG_INFINITY; k + 1]; n + 1];

    // Base cases
    for j in 0..=k {
        dp[0][j] = 0.0; // log(1) = 0, one way to partition 0
    }

    // Fill the dp table
    for i in 1..=n {
        for j in 1..=k {
            // Case 1: Don't use part j
            let log_without = dp[i][j - 1];

            // Case 2: Use at least one item in part j
            let log_with = if i >= j {
                dp[i - j][j]
            } else {
                f64::NEG_INFINITY // Not possible
            };

            // Combine using log-sum-exp
            if log_without == f64::NEG_INFINITY {
                dp[i][j] = log_with;
            } else if log_with == f64::NEG_INFINITY {
                dp[i][j] = log_without;
            } else {
                let max_log = log_without.max(log_with);
                let min_log = log_without.min(log_with);
                dp[i][j] = max_log + (min_log - max_log).exp().ln_1p();
            }
        }
    }

    dp[n][k]
}

/// Helper function to verify the algorithm's correctness for small inputs
fn count_configurations_exact(n: usize, k: usize, w: usize) -> Vec<usize> {
    let max_g = n / k;
    let mut results = Vec::with_capacity(max_g - w + 1);

    for g in w..=max_g {
        // Count partitions with exactly k parts, where the smallest part is g
        let count = count_partitions_with_min_part(n, k, g);
        results.push(count);
    }

    results
}

/// Recursively count partitions with specific constraints (for verification)
fn count_partitions_with_min_part(n: usize, k: usize, min_part: usize) -> usize {
    // Base cases
    if k == 1 {
        return if n >= min_part { 1 } else { 0 };
    }

    let mut count = 0;
    // Try different sizes for the largest part
    for first_part in min_part..=n {
        // Ensure we don't exceed remaining parts' capacity
        if n - first_part >= min_part * (k - 1) {
            count +=
                count_partitions_with_largest_part(n - first_part, k - 1, min_part, first_part);
        }
    }

    count
}

/// Count partitions where the largest part is exactly max_part
fn count_partitions_with_largest_part(
    n: usize,
    k: usize,
    min_part: usize,
    max_part: usize,
) -> usize {
    if k == 1 {
        return if n >= min_part && n <= max_part { 1 } else { 0 };
    }

    let mut count = 0;
    // Try different sizes for the next part, ensuring non-increasing order
    for next_part in min_part..=max_part.min(n) {
        if n - next_part >= min_part * (k - 1) {
            count += count_partitions_with_largest_part(n - next_part, k - 1, min_part, next_part);
        }
    }

    count
}

#[roxido]
fn count_for_size_configuration(n: usize, k: usize, w: usize) {
    count_configurations_log_iterative(n, k, w)
}

#[roxido]
fn count_for_size_configuration_old(n: usize, k: usize, w: usize) {
    let data = count_partitions_ln(n, k, w);
    let result = RMatrix::<f64>::new(2, data.len(), pc);
    let slice = result.slice_mut();
    for (i, d) in data.iter().enumerate() {
        slice[2 * i] = d.0 as f64;
        slice[2 * i + 1] = d.1;
    }
    result
}
