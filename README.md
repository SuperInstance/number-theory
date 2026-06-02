# numtheory

Number theory in Rust. Primes, residues, and the music of integers.

---

## What This Does

Pure Rust implementations of computational number theory:

- **Primes** — Sieve of Eratosthenes, trial division, deterministic Miller-Rabin (12 witnesses, correct for all u64), prime factorization, modular exponentiation, nth prime, prime counting function π(n).
- **Modular arithmetic** — Extended Euclidean algorithm, modular inverse, Chinese Remainder Theorem (CRT), Tonelli-Shanks modular square root, modular add/sub/mul.
- **Arithmetic functions** — Euler's totient φ(n) (single + sieve), Möbius function μ(n) (single + linear sieve), divisor count d(n), divisor sum σ(n) (overflow-safe u128 variant), Riemann zeta ζ(s) via Euler product, Mertens function M(n), gcd, lcm.
- **Continued fractions** — CF expansion of √n (periodic part extraction), convergent computation (h/k pairs), rational CF expansion, period length.
- **Quadratic residues** — Legendre symbol (a/p), Jacobi symbol (a/n), Kronecker symbol (a/n), quadratic residue testing, enumeration of all QRs mod p.
- **Diophantine equations** — Linear Diophantine solver (particular + general parameterized solution), Pell's equation x²−Dy²=1 (fundamental solution via continued fractions), negative Pell x²−Dy²=−1, nth Pell solution via exponentiation in ℤ[√D].
- **Dirichlet characters & L-functions** — Principal and Legendre characters, L(s,χ) evaluation for real and complex s, custom Complex type with arithmetic ops.

No external number theory dependencies — everything is implemented from scratch using `nalgebra`, `serde`, and the standard library.

---

## Install

```toml
[dependencies]
numtheory = "0.1"
```

Or clone and build:

```bash
git clone https://github.com/SuperInstance/number-theory.git
cd number-theory
cargo build
```

**Requirements:** Rust 2021 edition (≥ 1.56).

---

## Quick Start

### Prime generation and factorization

```rust
use numtheory::*;

let primes = sieve_primes(100);
println!("primes up to 100: {:?}", primes);

let factors = factorize(360);
for f in &factors {
    println!("{}^{}", f.prime, f.exponent);
}
// 2^3, 3^2, 5^1

assert!(is_prime(999999999999999989));
assert!(!is_prime(999999999999999990));
```

### Chinese Remainder Theorem

```rust
use numtheory::*;

// x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) → x = 23
let x = crt(&[2, 3, 2], &[3, 5, 7]);
assert_eq!(x, Some(23));
```

### Arithmetic functions

```rust
use numtheory::*;

println!("φ(30) = {}", euler_totient(30));  // 8
println!("μ(30) = {}", mobius(30));         // -1
println!("σ(28) = {}", divisor_sum_safe(28)); // 56 (perfect number)
println!("ζ(2) ≈ {:.6}", zeta_approx(2.0, 100)); // ≈ 1.6449... = π²/6
println!("M(100) = {}", mertens(100));       // 1
```

### Continued fractions and Pell's equation

```rust
use numtheory::*;

let cf = ContinuedFraction::from_sqrt(2);
println!("√2 = [{}; {:?}]", cf.a0, cf.periodic); // [1; [2]]

// Solve x² - 13y² = 1 → x=649, y=180
let sol = pell_solve(13).unwrap();
assert_eq!(sol.x, 649);
assert_eq!(sol.y, 180);
```

### Quadratic residues

```rust
use numtheory::*;

println!("(2/7) = {}", legendre(2, 7));   // 1 (QR)
println!("(3/7) = {}", legendre(3, 7));   // -1 (non-residue)
println!("QRs mod 7: {:?}", quadratic_residues(7)); // [1, 2, 4]
```

---

## API Reference

### Primes (`primes`)

| Function | Description |
|----------|-------------|
| `sieve_primes(limit)` | All primes ≤ limit via Eratosthenes |
| `is_prime_small(n)` | Trial division (n < 10⁶) |
| `miller_rabin(n, k)` | Miller-Rabin with 12 deterministic witnesses |
| `is_prime(n)` | Automatic: trial division or Miller-Rabin |
| `factorize(n)` | Prime factorization → `Vec<PrimeFactor>` |
| `mod_pow(base, exp, m)` | Modular exponentiation (u128 overflow-safe) |
| `nth_prime(n)` | 0-indexed nth prime |
| `prime_count(n)` | Prime counting function π(n) |

### Modular Arithmetic (`modular`)

| Function | Description |
|----------|-------------|
| `extended_gcd(a, b)` | Returns (gcd, x, y) with ax + by = gcd |
| `mod_inverse(a, m)` | Modular inverse via extended GCD |
| `crt(remainders, moduli)` | Chinese Remainder Theorem solver |
| `mod_sqrt(n, p)` | Tonelli-Shanks modular square root |
| `mod_add(a, b, m)`, `mod_sub(a, b, m)`, `mod_mul(a, b, m)` | Modular arithmetic |

### Arithmetic Functions (`arithmetic`)

| Function | Description |
|----------|-------------|
| `euler_totient(n)` | φ(n) — count of coprimes |
| `euler_totient_sieve(n)` | φ(1..=n) via sieve |
| `mobius(n)` | μ(n) — +1/−1/0 for square-free |
| `mobius_sieve(n)` | μ(1..=n) via linear sieve |
| `divisor_count(n)` | d(n) — number of divisors |
| `divisor_sum(n)` | σ(n) — sum of divisors (u64) |
| `divisor_sum_safe(n)` | σ(n) — overflow-safe (u128) |
| `zeta_approx(s, terms)` | ζ(s) via Euler product |
| `mertens(n)` | M(n) = Σ μ(k) for k=1..n |
| `gcd(a, b)`, `lcm(a, b)` | Greatest common divisor, least common multiple |

### Continued Fractions (`continued_fraction`)

| Method | Description |
|--------|-------------|
| `ContinuedFraction::from_sqrt(n)` | CF expansion of √n |
| `rational_cf(p, q)` | CF expansion of p/q |
| `.convergents(count)` | Best rational approximations (h, k) pairs |
| `.period_length()` | Period of the CF expansion |

### Quadratic Residues (`quadratic`)

| Function | Description |
|----------|-------------|
| `legendre(a, p)` | Legendre symbol (a/p) for odd prime p |
| `jacobi(a, n)` | Jacobi symbol (a/n) for odd n |
| `kronecker(a, n)` | Kronecker symbol (generalization) |
| `is_quadratic_residue(a, p)` | Boolean QR test |
| `quadratic_residues(p)` | All QRs mod p, sorted |

### Diophantine Equations (`diophantine`)

| Function | Description |
|----------|-------------|
| `linear_diophantine(a, b, c)` | Solve ax + by = c → particular solution |
| `linear_diophantine_general(a, b, c)` | Parameterized general solution |
| `pell_solve(d)` | Fundamental solution of x² − Dy² = 1 |
| `pell_solve_negative(d)` | Solve x² − Dy² = −1 |
| `pell_nth_solution(d, &fundamental, n)` | nth solution via √D-ring exponentiation |

### Dirichlet Characters & L-Functions (`dirichlet`)

| Method | Description |
|--------|-------------|
| `DirichletCharacter::principal(q)` | Principal character mod q |
| `DirichletCharacter::legendre_character(p)` | Legendre symbol character mod p |
| `.eval(n)` / `.eval_i64(n)` | Evaluate χ(n) |
| `l_function_at_1(&chi, terms)` | L(1, χ) via partial sums |
| `l_function(&chi, s, terms)` | L(s, χ) for complex s |
| `l_function_real(&chi, s, terms)` | L(s, χ) for real s |

---

## Tests

```bash
cargo test
```

---

## License

MIT OR Apache-2.0
