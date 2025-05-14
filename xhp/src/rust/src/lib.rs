// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

mod xhp;

use ahash::AHashMap;
use rand::prelude::*;
use rand::Rng;
use rand_distr::{Beta, BetaError, Binomial};
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

#[roxido]
fn entropy_from_partition(partition: &RVector) {
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
fn entropy_from_cluster_sizes(cluster_sizes: &RVector) {
    let cluster_sizes = cluster_sizes.to_i32(pc);
    let mut n_items = 0;
    let sum = cluster_sizes.slice().iter().fold(0.0, |s, &x| {
        let size = usize::try_from(x).stop();
        n_items += size;
        s - (size as f64) * (size as f64).ln()
    });
    sum / (n_items as f64) + (n_items as f64).ln()
}

pub fn sample_beta_binomial<R: Rng + ?Sized>(
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
