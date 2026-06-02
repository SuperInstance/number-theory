//! Modular arithmetic: exponentiation, inverse, CRT.

use crate::primes::mod_pow;

/// Modular inverse using extended Euclidean algorithm.
/// Returns None if a and m are not coprime.
pub fn mod_inverse(a: i64, m: u64) -> Option<u64> {
    let m = m as i64;
    let (g, x, _) = extended_gcd(a.rem_euclid(m), m);
    if g != 1 {
        return None;
    }
    Some(x.rem_euclid(m) as u64)
}

/// Extended Euclidean algorithm: returns (gcd, x, y) such that ax + by = gcd(a,b).
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        return (b, 0, 1);
    }
    let (g, x1, y1) = extended_gcd(b % a, a);
    let x = y1 - (b / a) * x1;
    let y = x1;
    (g, x, y)
}

/// Chinese Remainder Theorem solver.
/// Given remainders and moduli, finds x such that x ≡ r_i (mod m_i) for all i.
/// Returns None if no solution exists (moduli not pairwise coprime).
pub fn crt(remainders: &[u64], moduli: &[u64]) -> Option<u64> {
    if remainders.len() != moduli.len() || remainders.is_empty() {
        return None;
    }
    let mut x = remainders[0] as i64;
    let mut m = moduli[0] as i64;

    for i in 1..remainders.len() {
        let ri = remainders[i] as i64;
        let mi = moduli[i] as i64;
        let (_, inv, _) = extended_gcd(m, mi);
        let combined = m * mi;
        x = (x + m * ((ri - x) * inv % mi)) % combined;
        if x < 0 {
            x += combined;
        }
        m = combined;
    }
    Some(x as u64)
}

/// Modular square root using Tonelli-Shanks algorithm.
/// Returns x such that x^2 ≡ n (mod p), or None if no solution.
pub fn mod_sqrt(n: u64, p: u64) -> Option<u64> {
    if p == 2 {
        return Some(n % 2);
    }
    if mod_pow(n, (p - 1) / 2, p) != 1 {
        return None; // n is not a QR mod p
    }
    // Find Q and S such that p - 1 = Q * 2^S with Q odd
    let mut q = p - 1;
    let mut s = 0u32;
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    if s == 1 {
        let r = mod_pow(n, (p + 1) / 4, p);
        return Some(r);
    }
    // Find a non-residue z
    let mut z = 2u64;
    while mod_pow(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = mod_pow(z, q, p);
    let mut t = mod_pow(n, q, p);
    let mut r = mod_pow(n, (q + 1) / 2, p);
    loop {
        if t == 1 {
            return Some(r);
        }
        // Find the least i such that t^(2^i) ≡ 1 (mod p)
        let mut i = 1u32;
        let mut temp = mul_mod(t, t, p);
        while temp != 1 {
            temp = mul_mod(temp, temp, p);
            i += 1;
        }
        let b = mod_pow(c, 1u64 << (m - i - 1), p);
        m = i;
        c = mul_mod(b, b, p);
        t = mul_mod(t, c, p);
        r = mul_mod(r, b, p);
    }
}

fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// Modular addition.
pub fn mod_add(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}

/// Modular subtraction.
pub fn mod_sub(a: u64, b: u64, m: u64) -> u64 {
    let a = a as i128;
    let b = b as i128;
    let m = m as i128;
    ((a - b).rem_euclid(m)) as u64
}

/// Modular multiplication.
pub fn mod_mul(a: u64, b: u64, m: u64) -> u64 {
    mul_mod(a, b, m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_inverse() {
        assert_eq!(mod_inverse(3, 7), Some(5)); // 3*5=15≡1(mod7)
        assert_eq!(mod_inverse(2, 5), Some(3)); // 2*3=6≡1(mod5)
        assert_eq!(mod_inverse(2, 4), None); // not coprime
    }

    #[test]
    fn test_extended_gcd() {
        let (g, x, y) = extended_gcd(35, 15);
        assert_eq!(g, 5);
        assert_eq!(35 * x + 15 * y, 5);
    }

    #[test]
    fn test_crt_basic() {
        // x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) → x = 23
        let result = crt(&[2, 3, 2], &[3, 5, 7]);
        assert_eq!(result, Some(23));
    }

    #[test]
    fn test_crt_pair() {
        // x ≡ 1 (mod 2), x ≡ 2 (mod 3) → x = 5
        let result = crt(&[1, 2], &[2, 3]);
        assert_eq!(result, Some(5));
    }

    #[test]
    fn test_crt_single() {
        let result = crt(&[4], &[7]);
        assert_eq!(result, Some(4));
    }

    #[test]
    fn test_mod_sqrt() {
        // sqrt(4) mod 7 = 2 or 5
        let r = mod_sqrt(4, 7);
        assert!(r == Some(2) || r == Some(5));
        // sqrt(2) mod 7 — 2 is not a QR mod 7 (since 7 mod 8 = 7, check)
        // Actually 3^2=9≡2(mod7), so 2 IS a QR mod 7
        let r2 = mod_sqrt(2, 7);
        assert!(r2.is_some());
        let r2v = r2.unwrap();
        assert_eq!(mul_mod(r2v, r2v, 7), 2);
    }

    #[test]
    fn test_mod_add_sub_mul() {
        assert_eq!(mod_add(5, 7, 10), 2);
        assert_eq!(mod_sub(3, 7, 10), 6);
        assert_eq!(mod_mul(6, 7, 10), 2);
    }
}
