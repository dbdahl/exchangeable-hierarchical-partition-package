// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use roxido::*;
use std::f64;

#[roxido]
fn make_gupd(prob: &RVector, n_items: usize) {
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
            let dist = WeightedIndex::new(&weights).unwrap();
            let index = dist.sample(&mut rng);
            counts[index] += 1;
            r -= 1;
        }
        counts.sort();
        slice[i * n_clusters..(i + 1) * n_clusters].copy_from_slice(&counts);
    }
    result
}
