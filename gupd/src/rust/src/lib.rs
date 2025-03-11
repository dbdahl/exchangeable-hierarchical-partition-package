// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use ahash::AHashMap;
use rand::{Rng, SeedableRng};
use rand_distr::{Beta, BetaError, Binomial};
use rand_pcg::Pcg64Mcg;
use roxido::*;
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
    let sc = SizeConfiguration::sample(
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

    pub fn sample<R: Rng + ?Sized>(
        n_items: u64,
        n_clusters: u64,
        alpha: f64,
        beta: f64,
        rng: &mut R,
    ) -> Result<Self, &'static str> {
        let mut sg = Self::new_max_entropy(n_items, n_clusters)?;
        let mut index = sg.x.len() - 1;
        while index > 0 {
            let available = sg.available(index);
            let n = sample_beta_binomial(available, alpha, beta, rng)
                .map_err(|_| "Invalid beta parameter")?;
            sg.redistribute(index, n);
            index -= 1;
        }
        Ok(sg)
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
