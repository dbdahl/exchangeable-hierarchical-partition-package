// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();
use roxido::*;

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
    let sum = prob.iter().sum::<f64>();
    for p in prob.iter_mut() {
        *p /= sum;
    }
    let mut vec = vec![0; n_items + 1];
    vec[0] = 1; // Base case: S(0,0) = 1
    for i in 1..=n_items {
        for k in (1..=prob.len().min(i)).rev() {
            // Iterate backwards to avoid overwriting
            vec[k] = k * vec[k] + vec[k - 1]; // Recursive definition
        }
        vec[0] = 0; // Ensure S(n, 0) remains 0 for n > 0
    }
    let rval = RVector::from_value(0.0, prob.len(), pc);
    let slice = rval.slice_mut();
    for (r, (v, p)) in slice
        .iter_mut()
        .zip(vec.iter_mut().skip(1).zip(prob.iter()))
    {
        *r = *p / (*v as f64);
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
