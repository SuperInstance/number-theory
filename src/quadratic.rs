//! Quadratic residues and Legendre/Jacobi symbols.

use crate::primes::mod_pow;

/// Legendre symbol (a/p) for odd prime p.
/// Returns 1 if a is a QR mod p, -1 if not, 0 if p divides a.
pub fn legendre(a: i64, p: u64) -> i64 {
    let a = ((a % p as i64) + p as i64) as u64 % p;
    if a == 0 {
        return 0;
    }
    let result = mod_pow(a, (p - 1) / 2, p);
    if result == 1 {
        1
    } else {
        -1 // result == p-1
    }
}

/// Jacobi symbol (a/n) for odd positive integer n.
pub fn jacobi(mut a: i64, mut n: i64) -> i64 {
    assert!(n > 0 && n % 2 == 1, "n must be a positive odd integer");
    a = a % n;
    let mut result = 1i64;
    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            let n_mod_8 = n % 8;
            if n_mod_8 == 3 || n_mod_8 == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            result = -result;
        }
        a = a % n;
    }
    if n == 1 { result } else { 0 }
}

/// Kronecker symbol (a/n) — generalization of Jacobi symbol.
pub fn kronecker(a: i64, n: i64) -> i64 {
    if n == 0 {
        return if a == 1 || a == -1 { 1 } else { 0 };
    }

    let mut n = n;
    let a = a;

    // Handle n < 0
    let mut result = 1i64;
    if n < 0 {
        n = -n;
        if a < 0 {
            result = -1;
        }
    }

    // Handle even n
    if n % 2 == 0 {
        let v = n.trailing_zeros() as i64;
        n >>= v;
        if v % 2 == 1 && (a % 8 == 3 || a % 8 == 5) {
            result = -result;
        }
    }

    if n == 1 {
        return result;
    }

    result * jacobi(a, n)
}

/// Check if a is a quadratic residue mod p (p prime).
pub fn is_quadratic_residue(a: u64, p: u64) -> bool {
    legendre(a as i64, p) == 1
}

/// Find all quadratic residues mod p.
pub fn quadratic_residues(p: u64) -> Vec<u64> {
    let mut residues = Vec::new();
    for x in 1..p {
        let r = crate::primes::mod_pow(x, 2, p);
        if !residues.contains(&r) {
            residues.push(r);
        }
    }
    residues.sort();
    residues
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_legendre() {
        // (2/7) = 1 since 3^2 ≡ 2 (mod 7)
        assert_eq!(legendre(2, 7), 1);
        // (3/7) = -1 (non-residue... wait, let's check)
        // 3 is not a QR mod 7? 1^2=1, 2^2=4, 3^2=2, 4^2=2, 5^2=4, 6^2=1
        // QRs mod 7: {1, 2, 4}. So 3 is not a QR.
        // Wait, 5^2=25≡4, 6^2=36≡1. QRs = {1, 2, 4}. So 3 is not.
        assert_eq!(legendre(3, 7), -1);
        assert_eq!(legendre(7, 7), 0);
    }

    #[test]
    fn test_legendre_values() {
        // (5/11) — QRs mod 11: 1,3,4,5,9. So 5 is a QR.
        assert_eq!(legendre(5, 11), 1);
        // (2/11) — 2 is not a QR mod 11 (since 11 ≡ 3 mod 8)
        assert_eq!(legendre(2, 11), -1);
    }

    #[test]
    fn test_jacobi() {
        // (2/15) = (2/3)(2/5) = (-1)(-1) = 1
        assert_eq!(jacobi(2, 15), 1);
        // (5/21) = (5/3)(5/7) = (-1)(-1)... let's check
        // (5/3) = (2/3) = -1, (5/7) = (5/7) legendre = (5^3) mod 7 = 125 mod 7 = 6 = -1
        // So (5/21) = (-1)(-1) = 1
        assert_eq!(jacobi(5, 21), 1);
    }

    #[test]
    fn test_jacobi_prime_agrees() {
        // Jacobi with prime modulus should equal Legendre
        assert_eq!(jacobi(2, 7), legendre(2, 7));
        assert_eq!(jacobi(3, 7), legendre(3, 7));
    }

    #[test]
    fn test_is_quadratic_residue() {
        assert!(is_quadratic_residue(2, 7));
        assert!(!is_quadratic_residue(3, 7));
    }

    #[test]
    fn test_quadratic_residues() {
        let qrs = quadratic_residues(7);
        assert_eq!(qrs, vec![1, 2, 4]);
    }

    #[test]
    fn test_kronecker() {
        assert_eq!(kronecker(2, 15), 1);
        assert_eq!(kronecker(1, 0), 1);
    }
}
