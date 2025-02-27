// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use ibig::{ops::Abs, IBig, UBig};
use num_bigint::BigInt;
use num_traits::cast::ToPrimitive;
use num_traits::{One, Zero};
use roxido::*;
use rug::{Assign, Float}; // Added import for ibig

#[roxido]
fn make_gupd(prob: &RVector, n_items: usize, log: bool, method: &str) {
    println!("Method: {}", method);
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
    let rval = RVector::from_value(0.0, prob.len(), pc);
    let slice = rval.slice_mut();
    match method {
        "rug" => {
            let mut vec = vec![Float::with_val(256, 0); n_items + 1];
            vec[0] = Float::with_val(256, 1); // Base case: S(0,0) = 1
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
        }
        "num_bigint" => {
            let mut vec = vec![BigInt::from(0); n_items + 1];
            vec[0] = BigInt::one(); // Base case: S(0,0) = 1
            for i in 1..=n_items {
                for k in (1..=prob.len().min(i)).rev() {
                    vec[k] = &vec[k] * k + &vec[k - 1];
                }
                vec[0] = BigInt::from(0); // Ensure S(n, 0) remains 0 for n > 0
            }
            for (r, (v, p)) in slice.iter_mut().zip(vec.iter().skip(1).zip(prob.iter())) {
                if v.is_zero() {
                    stop!("Encountered zero in logarithm calculation.");
                } else {
                    *r = if log {
                        p.ln() - big_int_ln_precise(v)
                    } else {
                        p / v.to_f64().unwrap()
                    }
                }
            }
        }
        "ibig" => {
            let mut vec = vec![IBig::from(0); n_items + 1];
            vec[0] = IBig::from(1); // Base case: S(0,0) = 1
            for i in 1..=n_items {
                for k in (1..=prob.len().min(i)).rev() {
                    vec[k] = &vec[k] * k + &vec[k - 1];
                }
                vec[0] = IBig::from(0); // Ensure S(n, 0) remains 0 for n > 0
            }
            for (r, (v, p)) in slice.iter_mut().zip(vec.iter().skip(1).zip(prob.iter())) {
                if *v == IBig::from(0) {
                    stop!("Encountered zero in logarithm calculation.");
                } else {
                    *r = if log {
                        p.ln() - ibig_ln_precise(v)
                    } else {
                        p / ibig_to_f64(v)
                    }
                }
            }
        }
        _ => stop!("Unsupported method."),
    }
    rval
}

// Helper function to convert IBig to f64
fn ibig_to_f64(value: &IBig) -> f64 {
    // Try to convert to i64 first for small numbers
    if let Ok(i) = i64::try_from(value) {
        return i as f64;
    }

    // For large numbers, use the string representation
    let s = value.to_string();
    let sign = if s.starts_with('-') { -1.0 } else { 1.0 };

    // For very large numbers, approximate
    if s.len() > 18 {
        let digits = s.trim_start_matches('-').len() as f64;
        return sign * 10.0_f64.powf(digits - 1.0);
    }

    // Otherwise parse the string
    s.parse::<f64>().unwrap_or(f64::MAX * sign)
}

// Helper function to calculate natural logarithm of IBig with precision
fn ibig_ln_precise(value: &IBig) -> f64 {
    // Natural log approximation for large integers
    let s = value.abs().to_string();
    let length = s.len() as f64;

    if length <= 15.0 {
        // For smaller numbers that fit in f64
        return ibig_to_f64(value).ln();
    }

    // For large integers: ln(n) ≈ ln(10) * (number of digits - 1) + ln(leading digits)
    let leading_digits = s
        .chars()
        .take(15)
        .collect::<String>()
        .parse::<f64>()
        .unwrap();

    let ln_10 = 2.302585092994046; // ln(10)
    leading_digits.ln() + ln_10 * (length - 1.0)
}

// Function to calculate ln of BigInt with precision
fn big_int_ln_precise(n: &BigInt) -> f64 {
    if n <= &BigInt::from(0) {
        panic!("Cannot take logarithm of zero or negative number");
    }

    // For small numbers, direct conversion is fine
    if n <= &BigInt::from(u64::MAX) {
        return n.to_f64().unwrap().ln();
    }

    // Calculate the bit length (floor(log2(n)) + 1)
    let bit_length = n.bits();

    // Get the most significant ~53 bits
    let shift = if bit_length > 53 { bit_length - 53 } else { 0 };
    let leading_bits = if shift > 0 {
        (n >> shift).to_f64().unwrap()
    } else {
        n.to_f64().unwrap()
    };

    // Combine: ln(n) = ln(leading_bits * 2^shift) = ln(leading_bits) + shift*ln(2)
    leading_bits.ln() + (shift as f64) * std::f64::consts::LN_2
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
