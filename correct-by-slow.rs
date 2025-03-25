#!/usr/bin/env rust-script
//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! [dependencies]
//! rand = "0.9"
//! ```

use rand::distr::weighted::WeightedIndex;
use rand::prelude::*;

/// Computes the unsigned Stirling number of the first kind recursively,
/// returning the result as an f64.
///
/// Base cases:
/// - stirling(0,0) = 1
/// - stirling(n,0) = 0 or stirling(0,k) = 0 for n>0 or k>0
/// - stirling(n,n) = 1
///
/// Recursive case:
///   S(n,k) = S(n-1,k-1) + (n-1) * S(n-1,k)
fn stirling(n: usize, k: usize) -> f64 {
    if n == 0 && k == 0 {
        return 1.0;
    }
    if n == 0 || k == 0 {
        return 0.0;
    }
    if n == k {
        return 1.0;
    }
    stirling(n - 1, k - 1) + (n - 1) as f64 * stirling(n - 1, k)
}

/// Computes factorial(n) as an f64.
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, i| acc * (i as f64))
}

/// Computes the binomial coefficient "n choose k" as an f64.
/// Here we compute it as factorial(n) / (factorial(k) * factorial(n - k)).
fn choose(n: usize, k: usize) -> f64 {
    factorial(n) / (factorial(k) * factorial(n - k))
}

/// Recursively samples a cluster size configuration for n items partitioned into k clusters
/// using the following rule:
///
/// - If k == 1, return [n]
/// - Otherwise, let the first cluster size s be chosen from 1 ..= n - k + 1 with weight:
///       weight(s) = choose(n-1, s-1) * factorial(s-1) * stirling(n-s, k-1)
///   Then, recursively sample the sizes for the remaining clusters from the remaining items.
fn sample_cluster_sizes(n: usize, k: usize) -> Vec<usize> {
    if k == 1 {
        return vec![n];
    }

    let max_s = n - k + 1; // each remaining cluster must get at least one item
    let mut weights: Vec<f64> = Vec::with_capacity(max_s);
    let mut possible_s: Vec<usize> = Vec::with_capacity(max_s);

    // For each candidate first cluster size s, compute the weight.
    for s in 1..=max_s {
        let weight = choose(n - 1, s - 1) * factorial(s - 1) * stirling(n - s, k - 1);
        weights.push(weight);
        possible_s.push(s);
    }

    // Sample s from the possible sizes using the computed weights.
    let mut rng = rand::rng();
    let dist = WeightedIndex::new(&weights).expect("Weights should be nonnegative");
    let s_sample = possible_s[dist.sample(&mut rng)];

    // Recursively sample the remaining clusters.
    let mut result = vec![s_sample];
    let mut remaining = sample_cluster_sizes(n - s_sample, k - 1);
    result.append(&mut remaining);
    result
}

fn main() {
    // Example usage:
    // Let's say we have n = 10 items and we want exactly k = 3 clusters.
    let n = 10;
    let k = 3;
    let num_runs = 10_000_000;

    // HashMap to store frequency counts of the resulting strings.
    let mut freq: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // Run the sampling 10,000 times.
    for _ in 0..num_runs {
        let mut config = sample_cluster_sizes(n, k);
        // Sort the configuration in descending order.
        config.sort_unstable_by(|a, b| b.cmp(a));
        // Convert the configuration to a single string.
        let bob = config
            .iter()
            .map(|&x| x.to_string())
            .collect::<Vec<String>>()
            .join("");
        *freq.entry(bob).or_insert(0) += 1;
    }

    let mut results = freq.iter().collect::<Vec<_>>();
    results.sort_unstable_by(|a, b| a.cmp(b));
    // Print relative frequencies.
    for (s, count) in results.iter() {
        println!("{}: {:.4}", s, **count as f64 / num_runs as f64);
    }
}
