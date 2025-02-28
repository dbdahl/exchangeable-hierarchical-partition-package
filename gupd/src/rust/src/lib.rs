// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use ibig::UBig; // Added import for ibig
use num_bigint::BigInt;
use num_traits::cast::ToPrimitive;
use num_traits::{One, Zero};
use roxido::*;
use rug::{Assign, Float};

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
            // Using a similar approach to the rug implementation for better accuracy
            // Create a vector to store UBig values (using unsigned for slightly better performance)
            let mut vec = vec![UBig::from(0u32); n_items + 1];
            vec[0] = UBig::from(1u32); // Base case: S(0,0) = 1

            // Fill the vector using the recurrence relation
            for i in 1..=n_items {
                for k in (1..=prob.len().min(i)).rev() {
                    // S(n,k) = k*S(n-1,k) + S(n-1,k-1)
                    let k_ubig = UBig::from(k);
                    let term1 = &vec[k] * &k_ubig;
                    vec[k] = term1 + &vec[k - 1];
                }
                vec[0] = UBig::from(0u32); // Ensure S(n, 0) remains 0 for n > 0
            }

            // Compute the final result
            for (r, (v, p)) in slice.iter_mut().zip(vec.iter().skip(1).zip(prob.iter())) {
                if *v == UBig::from(0u32) {
                    stop!("Encountered zero in logarithm calculation.");
                } else {
                    *r = if log {
                        // For log calculation, use high-precision approach
                        p.ln() - ln_ubig_precise(v)
                    } else {
                        // For direct calculation
                        p / ubig_to_f64(v)
                    }
                }
            }
        }
        _ => stop!("Unsupported method."),
    }
    rval
}

// Helper function to convert UBig to f64 with high precision
fn ubig_to_f64(value: &UBig) -> f64 {
    // For small numbers, try to convert via u64
    if let Ok(small_val) = u64::try_from(value) {
        return small_val as f64;
    }

    // For larger numbers, convert via string
    let s = value.to_string();

    // If it's not too large, parse directly
    if s.len() <= 17 {
        return s.parse::<f64>().unwrap();
    }

    // For very large numbers, use scientific notation approach
    let digits = s.len();
    let leading_part = &s[0..std::cmp::min(16, s.len())];
    let leading_val: f64 = leading_part.parse().unwrap();

    // Compute as: leading_digits * 10^(remaining_digits)
    leading_val * 10f64.powf((digits - leading_part.len()) as f64)
}

// Helper function to calculate natural logarithm of UBig with high precision
fn ln_ubig_precise(value: &UBig) -> f64 {
    let s = value.to_string();

    // For small numbers, convert directly
    if s.len() <= 16 {
        return ubig_to_f64(value).ln();
    }

    // For large numbers, use ln(a * 10^b) = ln(a) + b * ln(10)
    let ln_10 = 2.302585092994046; // ln(10)

    // Get the most significant digits (as many as will fit in f64 without losing precision)
    let significant_digits = &s[0..std::cmp::min(15, s.len())];
    let significant_value: f64 = significant_digits.parse().unwrap();

    // Calculate ln(significant_value) + (total_digits - significant_digits) * ln(10)
    significant_value.ln() + ((s.len() - significant_digits.len()) as f64) * ln_10
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
