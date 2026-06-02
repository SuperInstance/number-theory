//! Continued fractions: convergents, quadratic irrationals.

use serde::{Deserialize, Serialize};

/// A continued fraction representation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContinuedFraction {
    /// Integer part (a0)
    pub a0: u64,
    /// Periodic part (empty for rationals)
    pub periodic: Vec<u64>,
}

impl ContinuedFraction {
    /// Compute continued fraction expansion of sqrt(n).
    pub fn from_sqrt(n: u64) -> Self {
        let sqrt_n = (n as f64).sqrt() as u64;
        if sqrt_n * sqrt_n == n {
            return ContinuedFraction { a0: sqrt_n, periodic: vec![] };
        }

        let a0 = sqrt_n;
        let mut periodic = Vec::new();
        let mut m = 0u64;
        let mut d = 1u64;
        let mut a = a0;

        loop {
            m = d * a - m;
            d = (n - m * m) / d;
            a = (a0 + m) / d;
            periodic.push(a);
            if a == 2 * a0 {
                break;
            }
        }
        ContinuedFraction { a0: a0, periodic }
    }

    /// Compute convergents as (numerator, denominator) pairs.
    /// Returns `count` convergents (or all for finite CF).
    pub fn convergents(&self, count: usize) -> Vec<(u64, u64)> {
        let mut terms = vec![self.a0];
        let mut i = 0;
        while terms.len() < count {
            terms.push(self.periodic[i % self.periodic.len()]);
            i += 1;
        }

        let mut convergents = Vec::new();
        let mut h_prev2 = 0u64;
        let mut h_prev1 = 1u64;
        let mut k_prev2 = 1u64;
        let mut k_prev1 = 0u64;

        for &a_i in &terms[..count.min(terms.len())] {
            let h = a_i * h_prev1 + h_prev2;
            let k = a_i * k_prev1 + k_prev2;
            convergents.push((h, k));
            h_prev2 = h_prev1;
            h_prev1 = h;
            k_prev2 = k_prev1;
            k_prev1 = k;
        }
        convergents
    }

    /// Get the period length.
    pub fn period_length(&self) -> usize {
        self.periodic.len()
    }
}

/// Compute continued fraction for a rational number p/q.
pub fn rational_cf(p: u64, q: u64) -> ContinuedFraction {
    if q == 0 {
        return ContinuedFraction { a0: p, periodic: vec![] };
    }
    let a0 = p / q;
    let mut terms = Vec::new();
    let mut rem = p % q;
    let mut d = q;
    while rem != 0 {
        let new_d = rem;
        let a = d / rem;
        rem = d % rem;
        d = new_d;
        terms.push(a);
    }
    ContinuedFraction { a0, periodic: terms }
}

/// Evaluate a continued fraction convergent as a float.
pub fn convergent_value(h: u64, k: u64) -> f64 {
    h as f64 / k as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cf_sqrt2() {
        let cf = ContinuedFraction::from_sqrt(2);
        assert_eq!(cf.a0, 1);
        assert_eq!(cf.periodic, vec![2]);
    }

    #[test]
    fn test_cf_sqrt3() {
        let cf = ContinuedFraction::from_sqrt(3);
        assert_eq!(cf.a0, 1);
        assert_eq!(cf.periodic, vec![1, 2]);
    }

    #[test]
    fn test_cf_sqrt7() {
        let cf = ContinuedFraction::from_sqrt(7);
        assert_eq!(cf.a0, 2);
        assert_eq!(cf.periodic, vec![1, 1, 1, 4]);
    }

    #[test]
    fn test_cf_perfect_square() {
        let cf = ContinuedFraction::from_sqrt(9);
        assert_eq!(cf.a0, 3);
        assert!(cf.periodic.is_empty());
    }

    #[test]
    fn test_convergents_sqrt2() {
        let cf = ContinuedFraction::from_sqrt(2);
        let convs = cf.convergents(6);
        // sqrt(2) ≈ 1.4142135...
        assert_eq!(convs[0], (1, 1));   // 1
        assert_eq!(convs[1], (3, 2));   // 1.5
        assert_eq!(convs[2], (7, 5));   // 1.4
        assert_eq!(convs[3], (17, 12)); // 1.4166...
        assert_eq!(convs[4], (41, 29)); // 1.4137...
        assert_eq!(convs[5], (99, 70)); // 1.41428...
    }

    #[test]
    fn test_convergents_approximate_sqrt2() {
        let cf = ContinuedFraction::from_sqrt(2);
        let convs = cf.convergents(15);
        let (h, k) = convs[14]; // 15th convergent
        let val = convergent_value(h, k);
        assert!((val - std::f64::consts::SQRT_2).abs() < 1e-10);
    }

    #[test]
    fn test_rational_cf() {
        let cf = rational_cf(22, 7);
        assert_eq!(cf.a0, 3);
        assert_eq!(cf.periodic, vec![7]);
    }

    #[test]
    fn test_rational_cf_simple() {
        let cf = rational_cf(7, 3);
        assert_eq!(cf.a0, 2);
        assert_eq!(cf.periodic, vec![3]);
    }

    #[test]
    fn test_period_length() {
        let cf = ContinuedFraction::from_sqrt(13);
        assert_eq!(cf.period_length(), 5);
    }
}
