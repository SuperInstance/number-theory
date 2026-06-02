//! Arithmetic functions: Euler's totient, Möbius, divisor sum, Riemann zeta approximation.

use crate::primes::{factorize, mod_pow, sieve_primes};

/// Euler's totient function φ(n): count of integers 1..n coprime to n.
pub fn euler_totient(n: u64) -> u64 {
    if n <= 1 {
        return n;
    }
    let factors = factorize(n);
    let mut result = n;
    for pf in &factors {
        result = result / pf.prime * (pf.prime - 1);
    }
    result
}

/// Euler's totient for all numbers 1..=n using sieve.
pub fn euler_totient_sieve(n: usize) -> Vec<u64> {
    let mut phi = (0..=n as u64).collect::<Vec<u64>>();
    for i in 2..=n {
        if phi[i] == i as u64 {
            // i is prime
            for j in (i..=n).step_by(i) {
                phi[j] = phi[j] / i as u64 * (i as u64 - 1);
            }
        }
    }
    phi
}

/// Möbius function μ(n):
/// - 1 if n is square-free with even number of prime factors
/// - -1 if n is square-free with odd number of prime factors
/// - 0 if n has a squared prime factor
pub fn mobius(n: u64) -> i64 {
    if n == 1 {
        return 1;
    }
    let factors = factorize(n);
    for pf in &factors {
        if pf.exponent > 1 {
            return 0;
        }
    }
    if factors.len() % 2 == 0 {
        1
    } else {
        -1
    }
}

/// Möbius function for all numbers 1..=n using sieve.
pub fn mobius_sieve(n: usize) -> Vec<i64> {
    let mut mu = vec![0i64; n + 1];
    mu[1] = 1;
    for i in 1..=n {
        if mu[i] == 0 {
            continue;
        }
        for j in (2 * i..=n).step_by(i) {
            mu[j] -= mu[i];
        }
    }
    // Zero out non-square-free numbers (they get wrong values from above)
    // Actually, the above sieve should be correct. Let's verify with factorize for small n.
    // The standard linear sieve is more correct. Let's use a simpler approach:
    let mut mu2 = vec![0i64; n + 1];
    mu2[1] = 1;
    let mut is_prime = vec![true; n + 1];
    let mut primes = Vec::new();
    for i in 2..=n {
        if is_prime[i] {
            primes.push(i);
            mu2[i] = -1;
        }
        for &p in &primes {
            if i * p > n {
                break;
            }
            is_prime[i * p] = false;
            if i % p == 0 {
                mu2[i * p] = 0;
                break;
            } else {
                mu2[i * p] = -mu2[i];
            }
        }
    }
    mu2
}

/// Number of divisors d(n).
pub fn divisor_count(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let factors = factorize(n);
    factors.iter().map(|pf| pf.exponent as u64 + 1).product()
}

/// Sum of divisors σ(n).
pub fn divisor_sum(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let factors = factorize(n);
    factors.iter().map(|pf| {
        let p = pf.prime;
        let k = pf.exponent as u64;
        // (p^(k+1) - 1) / (p - 1)
        (mod_pow(p, k + 1, u64::MAX).saturating_sub(1)) / (p - 1)
    }).product()
}

/// Sum of divisors with overflow-safe computation.
pub fn divisor_sum_safe(n: u64) -> u128 {
    if n == 0 {
        return 0;
    }
    let factors = factorize(n);
    factors.iter().map(|pf| {
        let p = pf.prime as u128;
        let k = pf.exponent as u32;
        (p.pow(k + 1) - 1) / (p - 1)
    }).product()
}

/// Riemann zeta function approximation ζ(s) using Euler product.
/// Uses the first `terms` primes.
pub fn zeta_approx(s: f64, terms: usize) -> f64 {
    let primes = sieve_primes(terms * 20); // enough primes
    let primes = &primes[..terms.min(primes.len())];
    let mut product = 1.0;
    for &p in primes {
        let ps = (p as f64).powf(s);
        product *= ps / (ps - 1.0);
    }
    product
}

/// Mertens function M(n) = sum of μ(k) for k = 1..n.
pub fn mertens(n: usize) -> i64 {
    let mu = mobius_sieve(n);
    mu[1..=n].iter().sum()
}

/// GCD of two numbers.
pub fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// LCM of two numbers.
pub fn lcm(a: u64, b: u64) -> u64 {
    a / gcd(a, b) * b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_totient() {
        assert_eq!(euler_totient(1), 1);
        assert_eq!(euler_totient(7), 6); // prime
        assert_eq!(euler_totient(10), 4); // 1,3,7,9
        assert_eq!(euler_totient(12), 4); // 1,5,7,11
        assert_eq!(euler_totient(30), 8);
    }

    #[test]
    fn test_euler_totient_sieve() {
        let phi = euler_totient_sieve(30);
        assert_eq!(phi[1], 1);
        assert_eq!(phi[7], 6);
        assert_eq!(phi[10], 4);
        assert_eq!(phi[12], 4);
        assert_eq!(phi[30], 8);
    }

    #[test]
    fn test_mobius() {
        assert_eq!(mobius(1), 1);
        assert_eq!(mobius(2), -1); // one prime factor
        assert_eq!(mobius(6), 1); // 2*3, two primes
        assert_eq!(mobius(4), 0); // 2^2, squared factor
        assert_eq!(mobius(30), -1); // 2*3*5, three primes
        assert_eq!(mobius(12), 0); // 2^2*3
    }

    #[test]
    fn test_mobius_sieve() {
        let mu = mobius_sieve(30);
        assert_eq!(mu[1], 1);
        assert_eq!(mu[2], -1);
        assert_eq!(mu[6], 1);
        assert_eq!(mu[4], 0);
        assert_eq!(mu[30], -1);
    }

    #[test]
    fn test_divisor_count() {
        assert_eq!(divisor_count(1), 1);
        assert_eq!(divisor_count(6), 4); // 1,2,3,6
        assert_eq!(divisor_count(12), 6); // 1,2,3,4,6,12
        assert_eq!(divisor_count(28), 6); // perfect number
    }

    #[test]
    fn test_divisor_sum() {
        assert_eq!(divisor_sum_safe(6), 12); // 1+2+3+6=12
        assert_eq!(divisor_sum_safe(28), 56); // 1+2+4+7+14+28=56
        assert_eq!(divisor_sum_safe(12), 28); // 1+2+3+4+6+12=28
    }

    #[test]
    fn test_zeta_approx() {
        let z2 = zeta_approx(2.0, 100);
        assert!((z2 - (std::f64::consts::PI * std::f64::consts::PI / 6.0)).abs() < 0.01);
    }

    #[test]
    fn test_gcd_lcm() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(lcm(12, 8), 24);
        assert_eq!(gcd(17, 13), 1);
    }
}
