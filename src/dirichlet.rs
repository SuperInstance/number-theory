//! Dirichlet characters and L-functions (basic).

use serde::{Deserialize, Serialize};

/// A Dirichlet character mod q.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DirichletCharacter {
    modulus: u64,
    /// Mapping from residue class to character value (as complex-like i64 pair: (re*SCALE, im*SCALE))
    values: Vec<i64>,
}

/// Complex value for character values (simplified as f64 pairs).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }

    pub fn one() -> Self {
        Complex { re: 1.0, im: 0.0 }
    }

    pub fn zero() -> Self {
        Complex { re: 0.0, im: 0.0 }
    }

    pub fn norm(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    pub fn pow(self, n: u64) -> Complex {
        if n == 0 {
            return Complex::one();
        }
        let mut result = Complex::one();
        let mut base = self;
        let mut exp = n;
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base;
            }
            base = base * base;
            exp /= 2;
        }
        result
    }

    pub fn exp(z: Complex) -> Complex {
        let e_re = z.re.exp();
        Complex {
            re: e_re * z.im.cos(),
            im: e_re * z.im.sin(),
        }
    }

    pub fn ln(&self) -> Complex {
        Complex {
            re: self.norm().ln(),
            im: self.re.atan2(self.im),
        }
    }
}

impl std::ops::Add for Complex {
    type Output = Complex;
    fn add(self, other: Complex) -> Complex {
        Complex { re: self.re + other.re, im: self.im + other.im }
    }
}

impl std::ops::Sub for Complex {
    type Output = Complex;
    fn sub(self, other: Complex) -> Complex {
        Complex { re: self.re - other.re, im: self.im - other.im }
    }
}

impl std::ops::Mul for Complex {
    type Output = Complex;
    fn mul(self, other: Complex) -> Complex {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl std::ops::Div for Complex {
    type Output = Complex;
    fn div(self, other: Complex) -> Complex {
        let denom = other.re * other.re + other.im * other.im;
        Complex {
            re: (self.re * other.re + self.im * other.im) / denom,
            im: (self.im * other.re - self.re * other.im) / denom,
        }
    }
}

impl DirichletCharacter {
    /// Create the principal character mod q (χ(n) = 1 if gcd(n,q)=1, else 0).
    pub fn principal(modulus: u64) -> Self {
        let mut values = vec![0i64; modulus as usize];
        for i in 1..modulus as usize {
            if crate::arithmetic::gcd(i as u64, modulus) == 1 {
                values[i] = 1;
            }
        }
        DirichletCharacter { modulus, values }
    }

    /// Create the Legendre symbol character mod p (for odd prime p).
    pub fn legendre_character(p: u64) -> Self {
        let mut values = vec![0i64; p as usize];
        for n in 1..p as usize {
            values[n] = crate::quadratic::legendre(n as i64, p);
        }
        DirichletCharacter { modulus: p, values }
    }

    /// Evaluate the character at n.
    pub fn eval_i64(&self, n: u64) -> i64 {
        if n == 0 {
            return 0;
        }
        self.values[(n % self.modulus) as usize]
    }

    /// Evaluate as complex number.
    pub fn eval(&self, n: u64) -> Complex {
        let v = self.eval_i64(n);
        Complex::new(v as f64, 0.0)
    }

    /// The modulus of this character.
    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Check if this is the principal character.
    pub fn is_principal(&self) -> bool {
        self.values.iter().filter(|&&v| v != 0).all(|&v| v == 1)
    }
}

/// Compute L(1, χ) using partial sums.
/// L(s, χ) = sum_{n=1}^{∞} χ(n) / n^s
pub fn l_function_at_1(chi: &DirichletCharacter, terms: usize) -> f64 {
    let mut sum = 0.0;
    for n in 1..=terms {
        let chi_n = chi.eval_i64(n as u64) as f64;
        sum += chi_n / n as f64;
    }
    sum
}

/// Compute L(s, χ) for complex s using partial sums.
pub fn l_function(chi: &DirichletCharacter, s: Complex, terms: usize) -> Complex {
    let mut sum = Complex::zero();
    for n in 1..=terms {
        let chi_n = chi.eval(n as u64);
        let n_complex = Complex::new(n as f64, 0.0);
        let n_s = Complex::exp(s * n_complex.ln());
        sum = sum + chi_n / n_s;
    }
    sum
}

/// Compute Dirichlet L-function approximation for real s.
pub fn l_function_real(chi: &DirichletCharacter, s: f64, terms: usize) -> f64 {
    let mut sum = 0.0;
    for n in 1..=terms {
        let chi_n = chi.eval_i64(n as u64) as f64;
        sum += chi_n / (n as f64).powf(s);
    }
    sum
}

/// Count primitive Dirichlet characters mod q.
pub fn count_primitive_characters(q: u64) -> u64 {
    crate::arithmetic::euler_totient(q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_principal_character() {
        let chi = DirichletCharacter::principal(5);
        assert_eq!(chi.eval_i64(1), 1);
        assert_eq!(chi.eval_i64(2), 1);
        assert_eq!(chi.eval_i64(3), 1);
        assert_eq!(chi.eval_i64(4), 1);
        assert_eq!(chi.eval_i64(5), 0); // gcd(5,5)≠1
        assert!(chi.is_principal());
    }

    #[test]
    fn test_legendre_character_mod5() {
        let chi = DirichletCharacter::legendre_character(5);
        // QRs mod 5: 1, 4. So (1/5)=1, (2/5)=-1, (3/5)=-1, (4/5)=1
        assert_eq!(chi.eval_i64(1), 1);
        assert_eq!(chi.eval_i64(2), -1);
        assert_eq!(chi.eval_i64(3), -1);
        assert_eq!(chi.eval_i64(4), 1);
        assert!(!chi.is_principal());
    }

    #[test]
    fn test_l_function_principal() {
        // L(s, χ₀) for principal character mod q approaches
        // φ(q)/q * ζ(s) as terms → ∞
        let chi = DirichletCharacter::principal(3);
        let l1 = l_function_at_1(&chi, 10000);
        // Should be close to φ(3)/3 * ∞... actually L(1, principal mod q) diverges
        // like harmonic series scaled. Let's just check it's positive and reasonable.
        assert!(l1 > 0.5);
    }

    #[test]
    fn test_l_function_legendre() {
        let chi = DirichletCharacter::legendre_character(5);
        let l1 = l_function_at_1(&chi, 10000);
        // L(1, legendre mod 5) should converge to something nonzero
        assert!(l1.abs() > 0.1);
        assert!(l1.abs() < 5.0);
    }

    #[test]
    fn test_l_function_real() {
        let chi = DirichletCharacter::principal(3);
        let l2 = l_function_real(&chi, 2.0, 1000);
        // Should be positive and finite
        assert!(l2 > 0.0);
    }

    #[test]
    fn test_complex_ops() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        let c = a * b;
        // (1+2i)(3+4i) = 3+4i+6i+8i² = 3+10i-8 = -5+10i
        assert!((c.re - (-5.0)).abs() < 1e-10);
        assert!((c.im - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_complex_pow() {
        let i = Complex::new(0.0, 1.0);
        let i2 = i.pow(2);
        assert!((i2.re - (-1.0)).abs() < 1e-10);
        assert!(i2.im.abs() < 1e-10);
    }
}
