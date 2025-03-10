// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use ahash::AHashMap;
use rand_distr::{Beta, BetaError, Binomial};
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
fn initial_state_r(n_items: i32, n_clusters: i32) {
    let r1 = SizeConfiguration::new_max_entropy(
        u64::try_from(n_items).stop(),
        u64::try_from(n_clusters).stop(),
    );
    let r2: Vec<_> = r1
        .stop()
        .x
        .iter()
        .map(|&x| i32::try_from(x).stop())
        .collect();
    r2
}

#[roxido]
fn new_from_vec(x: &RVector) {
    let x = x.to_i32(pc);
    let y = x
        .slice()
        .iter()
        .map(|&yy| u64::try_from(yy - 1).stop_str("Cluster sizes must be at least one"))
        .collect();
    let z = SizeConfiguration::new_from_vec(y).stop();
    RExternalPtr::encode(z, "SizeConfiguration", pc)
}

#[roxido]
fn new_max_entropy(n_items: i32, n_clusters: i32) {
    let y = SizeConfiguration::new_max_entropy(
        u64::try_from(n_items).stop(),
        u64::try_from(n_clusters).stop(),
    )
    .stop();
    RExternalPtr::encode(y, "SizeConfiguration", pc)
}

#[roxido]
fn size_configuration_to_r(x: &RExternalPtr) {
    let r1 = x.decode_ref::<SizeConfiguration>();
    let r2: Vec<_> = r1.x.iter().map(|&x| i32::try_from(x).stop()).collect();
    r2
}

#[roxido]
fn size_configuration_available(x: &RExternalPtr, index: i32) {
    let x = x.decode_ref::<SizeConfiguration>();
    i32::try_from(x.available(usize::try_from(index - 1).stop()).n).stop()
}

#[roxido]
fn size_configuration_redistribute(x: &mut RExternalPtr, index: i32, n: i32) {
    let x = x.decode_mut::<SizeConfiguration>();
    x.redistribute(
        usize::try_from(index - 1).stop(),
        Available {
            n: u64::try_from(n).stop(),
        },
    )
}

#[derive(Debug)]
pub struct Available {
    n: u64,
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
    fn new_from_vec(mut x: Vec<u64>) -> Result<Self, &'static str> {
        x.sort_unstable_by(|a, b| b.cmp(a));
        Ok(Self { x })
    }

    fn new_max_entropy(n_items: u64, n_clusters: u64) -> Result<Self, &'static str> {
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

    fn available(&self, index: usize) -> Available {
        Available {
            n: match index {
                0 => 0,
                i if i >= self.x.len() => 0,
                i if i == self.x.len() - 1 => self.x[i],
                i => self.x[i] - self.x[i + 1],
            },
        }
    }

    fn redistribute(&mut self, index: usize, mut n: Available) -> bool {
        if n.n == 0 {
            return true;
        }
        if index == 0 {
            return false;
        }
        self.x[index] -= n.n;
        let v0 = self.x[0];
        let mut i = index - 1;
        while self.x[i] != v0 {
            i -= 1;
        }
        i += 1;
        while i < index && n.n > 0 {
            self.x[i] += 1;
            i += 1;
            n.n -= 1;
        }
        if n.n == 0 {
            return true;
        }
        let whole = n.n / u64::try_from(index).expect("usize doesn't fit in u64");
        if whole > 0 {
            for j in 0..index {
                self.x[j] += whole;
            }
            n.n -= u64::try_from(index).expect("usize doesn't fit in u64") * whole;
        }
        for j in 0..usize::try_from(n.n).expect("u64 doesn't fit in usize") {
            self.x[j] += 1;
        }
        true
    }
}
