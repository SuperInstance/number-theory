//! Prime generation and primality testing.

use serde::{Deserialize, Serialize};

/// Sieve of Eratosthenes: generate all primes up to `limit`.
pub fn sieve_primes(limit: usize) -> Vec<u64> {
    if limit < 2 {
        return vec![];
    }
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;
    let mut i = 2;
    while i * i <= limit {
        if is_prime[i] {
            let mut j = i * i;
            while j <= limit {
                is_prime[j] = false;
                j += i;
            }
        }
        i += 1;
    }
    (2..=limit).filter(|&i| is_prime[i]).map(|i| i as u64).collect()
}

/// Deterministic primality test for small numbers using trial division.
pub fn is_prime_small(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n < 4 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    let mut i = 5u64;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

/// Miller-Rabin primality test (probabilistic).
/// Uses `k` rounds of testing. For deterministic results up to 3.3 × 10^24,
/// use the witnesses [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37].
pub fn miller_rabin(n: u64, _k: usize) -> bool {
    if n < 2 {
        return false;
    }
    if n < 4 {
        return true;
    }
    if n % 2 == 0 {
        return false;
    }

    // Write n-1 as 2^r * d
    let mut d = n - 1;
    let mut r = 0u32;
    while d % 2 == 0 {
        d /= 2;
        r += 1;
    }

    let witnesses = vec![2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    // These 12 witnesses are deterministic for all n < 3.3 × 10^24
    // Since u64 max < that threshold, this is always deterministic

    'witness: for &a in &witnesses {
        if a >= n {
            continue;
        }
        let mut x = mod_pow(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..r - 1 {
            x = mod_pow(x, 2, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// General primality test: uses deterministic Miller-Rabin for all practical sizes.
pub fn is_prime(n: u64) -> bool {
    if n < 1_000_000 {
        return is_prime_small(n);
    }
    miller_rabin(n, 12)
}

/// Factorize n into prime factors with their exponents.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimeFactor {
    pub prime: u64,
    pub exponent: u32,
}

/// Return the prime factorization of n.
pub fn factorize(mut n: u64) -> Vec<PrimeFactor> {
    if n <= 1 {
        return vec![];
    }
    let mut factors = Vec::new();
    let mut count = 0u32;
    while n % 2 == 0 {
        count += 1;
        n /= 2;
    }
    if count > 0 {
        factors.push(PrimeFactor { prime: 2, exponent: count });
    }
    let mut p = 3u64;
    while p * p <= n {
        let mut count = 0u32;
        while n % p == 0 {
            count += 1;
            n /= p;
        }
        if count > 0 {
            factors.push(PrimeFactor { prime: p, exponent: count });
        }
        p += 2;
    }
    if n > 1 {
        factors.push(PrimeFactor { prime: n, exponent: 1 });
    }
    factors
}

/// Modular exponentiation: computes base^exp mod m.
pub fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1u64;
    base %= modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = mul_mod(result, base, modulus);
        }
        exp /= 2;
        base = mul_mod(base, base, modulus);
    }
    result
}

/// Modular multiplication avoiding overflow.
fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Generate the nth prime (0-indexed). Uses sieve with growing bounds.
pub fn nth_prime(n: usize) -> u64 {
    if n == 0 {
        return 2;
    }
    // Upper bound approximation: p_n < n * (ln(n) + ln(ln(n))) for n >= 6
    let limit = if n < 6 {
        20
    } else {
        let nn = n as f64;
        (nn * (nn.ln() + nn.ln().ln()) * 1.3) as usize
    };
    let primes = sieve_primes(limit);
    primes[n]
}

/// Count primes up to n (prime counting function π(n)).
pub fn prime_count(n: u64) -> u64 {
    sieve_primes(n as usize).len() as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sieve_small() {
        let primes = sieve_primes(30);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn test_sieve_empty() {
        assert!(sieve_primes(1).is_empty());
    }

    #[test]
    fn test_is_prime_small() {
        assert!(!is_prime_small(0));
        assert!(!is_prime_small(1));
        assert!(is_prime_small(2));
        assert!(is_prime_small(3));
        assert!(!is_prime_small(4));
        assert!(is_prime_small(5));
        assert!(is_prime_small(97));
        assert!(!is_prime_small(99));
    }

    #[test]
    fn test_miller_rabin() {
        // Known primes
        assert!(miller_rabin(2, 5));
        assert!(miller_rabin(7919, 5));
        assert!(miller_rabin(104729, 5));
        assert!(miller_rabin(999999999999999989, 12));

        // Known composites
        assert!(!miller_rabin(4, 5));
        assert!(!miller_rabin(561, 12)); // Carmichael number
        assert!(!miller_rabin(100, 5));
    }

    #[test]
    fn test_is_prime_large() {
        assert!(is_prime(999999999999999989));
        assert!(!is_prime(999999999999999990));
    }

    #[test]
    fn test_factorize() {
        let f = factorize(60);
        assert_eq!(f.len(), 3);
        assert_eq!(f[0].prime, 2);
        assert_eq!(f[0].exponent, 2);
        assert_eq!(f[1].prime, 3);
        assert_eq!(f[1].exponent, 1);
        assert_eq!(f[2].prime, 5);
        assert_eq!(f[2].exponent, 1);
    }

    #[test]
    fn test_factorize_prime() {
        let f = factorize(17);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].prime, 17);
        assert_eq!(f[0].exponent, 1);
    }

    #[test]
    fn test_factorize_one() {
        assert!(factorize(1).is_empty());
    }

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 10, 1000), 1024 % 1000);
        assert_eq!(mod_pow(3, 7, 13), ((3u64.pow(7)) % 13));
        assert_eq!(mod_pow(5, 0, 7), 1);
    }

    #[test]
    fn test_nth_prime() {
        assert_eq!(nth_prime(0), 2);
        assert_eq!(nth_prime(1), 3);
        assert_eq!(nth_prime(4), 11);
        assert_eq!(nth_prime(9), 29);
    }

    #[test]
    fn test_prime_count() {
        assert_eq!(prime_count(10), 4);
        assert_eq!(prime_count(100), 25);
    }
}
