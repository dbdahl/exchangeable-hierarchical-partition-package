// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use roxido::*;
use rug::{Assign, Float};

#[roxido]
fn make_gupd(prob: &RVector, n_items: usize, log: bool) {
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

    // Use Float instead of BigUint
    let mut vec = vec![Float::with_val(256, 0); n_items + 1];
    vec[0].assign(Float::with_val(256, 1)); // Base case: S(0,0) = 1

    for i in 1..=n_items {
        for k in (1..=prob.len().min(i)).rev() {
            let k_float = Float::with_val(256, k);
            let mut temp = vec[k].clone();
            temp *= &k_float;
            temp += &vec[k - 1];
            vec[k] = temp;
        }
        vec[0].assign(Float::with_val(256, 0)); // Ensure S(n, 0) remains 0 for n > 0
    }

    let rval = RVector::from_value(0.0, prob.len(), pc);
    let slice = rval.slice_mut();

    for (r, (v, p)) in slice.iter_mut().zip(vec.iter().skip(1).zip(prob.iter())) {
        if v.is_zero() {
            stop!("Encountered zero in logarithm calculation.");
        } else {
            *r = if log {
                (Float::with_val(256, p.ln()) - v.clone().ln()).to_f64()
            } else {
                (Float::with_val(256, p) / v.clone()).to_f64()
            }
        }
    }

    rval
}

#[allow(dead_code)]
fn stirling_second_kind(n: usize, k: usize) -> usize {
    if n == 0 && k == 0 {
        1
    } else if k == 0 || n == 0 || k > n {
        0
    } else {
        k * stirling_second_kind(n - 1, k) + stirling_second_kind(n - 1, k - 1)
    }
}
