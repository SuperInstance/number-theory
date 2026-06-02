//! Diophantine equations: linear, Pell's equation.

use crate::modular::extended_gcd;
use crate::continued_fraction::ContinuedFraction;
use serde::{Deserialize, Serialize};

/// Solve linear Diophantine equation ax + by = c.
/// Returns Some((x, y)) if a solution exists, None otherwise.
pub fn linear_diophantine(a: i64, b: i64, c: i64) -> Option<(i64, i64)> {
    let (g, x, y) = extended_gcd(a.abs(), b.abs());
    if c % g != 0 {
        return None;
    }
    let factor = c / g;
    let x = x * factor * a.signum();
    let y = y * factor * b.signum();
    Some((x, y))
}

/// General solution of linear Diophantine equation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearDiophantineSolution {
    pub particular: (i64, i64), // (x0, y0)
    pub homogeneous: (i64, i64), // (dx, dy) — general: x = x0 + t*dx, y = y0 + t*dy
}

/// Solve with general solution parameterized by t ∈ Z.
pub fn linear_diophantine_general(a: i64, b: i64, c: i64) -> Option<LinearDiophantineSolution> {
    let (x0, y0) = linear_diophantine(a, b, c)?;
    let g = extended_gcd(a.abs(), b.abs()).0;
    let dx = b.abs() / g;
    let dy = -(a.abs() / g);
    Some(LinearDiophantineSolution {
        particular: (x0, y0),
        homogeneous: (dx, dy),
    })
}

/// Solution to Pell's equation x² - Dy² = 1.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PellSolution {
    pub x: u64,
    pub y: u64,
}

/// Solve Pell's equation x² - Dy² = 1 for non-square D > 1.
/// Returns the fundamental (minimal) solution using continued fractions.
pub fn pell_solve(d: u64) -> Option<PellSolution> {
    // Check D is not a perfect square
    let sqrt_d = (d as f64).sqrt() as u64;
    if sqrt_d * sqrt_d == d {
        return None;
    }
    if d <= 1 {
        return None;
    }

    let cf = ContinuedFraction::from_sqrt(d);
    let period = cf.period_length();

    // For x² - Dy² = 1, we need convergents at position:
    // - period-1 if period is odd (actually period if odd for x²-Dy²=1)
    // - 2*period-1 if period is even
    // Wait, let me be more careful:
    // The fundamental solution to x² - Dy² = 1 is the convergent at index:
    //   r-1 where r = period length, if r is even
    //   2r-1 if r is odd
    // Actually, the standard approach: try convergents until we find one that satisfies.

    let num_convergents = if period % 2 == 0 { period } else { 2 * period };
    let convs = cf.convergents(num_convergents);

    for (h, k) in &convs {
        // Check h² - D*k² = 1
        let lhs = *h as i128 * *h as i128 - d as i128 * *k as i128 * *k as i128;
        if lhs == 1 {
            return Some(PellSolution { x: *h, y: *k });
        }
    }

    // Fallback: try more convergents
    let convs = cf.convergents(4 * period);
    for (h, k) in &convs {
        let lhs = *h as i128 * *h as i128 - d as i128 * *k as i128 * *k as i128;
        if lhs == 1 {
            return Some(PellSolution { x: *h, y: *k });
        }
    }

    None
}

/// Solve the negative Pell's equation x² - Dy² = -1.
pub fn pell_solve_negative(d: u64) -> Option<PellSolution> {
    let sqrt_d = (d as f64).sqrt() as u64;
    if sqrt_d * sqrt_d == d {
        return None;
    }

    let cf = ContinuedFraction::from_sqrt(d);
    let period = cf.period_length();

    // x² - Dy² = -1 has a solution iff the period is odd
    if period % 2 == 0 {
        return None;
    }

    let convs = cf.convergents(2 * period);
    for (h, k) in &convs {
        let lhs = *h as i128 * *h as i128 - d as i128 * *k as i128 * *k as i128;
        if lhs == -1 {
            return Some(PellSolution { x: *h, y: *k });
        }
    }

    None
}

/// Generate further Pell solutions from fundamental one.
/// If (x1, y1) is the fundamental solution, the nth solution is:
/// x_n + y_n√D = (x1 + y1√D)^n
pub fn pell_nth_solution(d: u64, fundamental: &PellSolution, n: u32) -> PellSolution {
    let mut x = 1u128;
    let mut y = 0u128;
    let mut bx = fundamental.x as u128;
    let mut by = fundamental.y as u128;

    let mut exp = n;
    while exp > 0 {
        if exp % 2 == 1 {
            let new_x = x * bx + d as u128 * y * by;
            let new_y = x * by + y * bx;
            x = new_x;
            y = new_y;
        }
        let new_bx = bx * bx + d as u128 * by * by;
        let new_by = 2 * bx * by;
        bx = new_bx;
        by = new_by;
        exp /= 2;
    }

    PellSolution { x: x as u64, y: y as u64 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_diophantine() {
        // 3x + 5y = 8 → x=1, y=1 is a solution
        let sol = linear_diophantine(3, 5, 8);
        assert!(sol.is_some());
        let (x, y) = sol.unwrap();
        assert_eq!(3 * x + 5 * y, 8);
    }

    #[test]
    fn test_linear_diophantine_no_solution() {
        // 2x + 4y = 5 → no solution (gcd(2,4)=2 does not divide 5)
        assert!(linear_diophantine(2, 4, 5).is_none());
    }

    #[test]
    fn test_linear_diophantine_general() {
        let sol = linear_diophantine_general(6, 9, 21);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        let (x0, y0) = sol.particular;
        assert_eq!(6 * x0 + 9 * y0, 21);
    }

    #[test]
    fn test_pell_d2() {
        // x² - 2y² = 1 → fundamental: x=3, y=2
        let sol = pell_solve(2);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert_eq!(sol.x, 3);
        assert_eq!(sol.y, 2);
    }

    #[test]
    fn test_pell_d3() {
        // x² - 3y² = 1 → fundamental: x=2, y=1
        let sol = pell_solve(3);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert_eq!(sol.x, 2);
        assert_eq!(sol.y, 1);
    }

    #[test]
    fn test_pell_d5() {
        // x² - 5y² = 1 → fundamental: x=9, y=4
        let sol = pell_solve(5);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert_eq!(sol.x as u128 * sol.x as u128 - 5u128 * sol.y as u128 * sol.y as u128, 1);
        assert_eq!(sol.x, 9);
        assert_eq!(sol.y, 4);
    }

    #[test]
    fn test_pell_d7() {
        // x² - 7y² = 1 → fundamental: x=8, y=3
        let sol = pell_solve(7);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert_eq!(sol.x, 8);
        assert_eq!(sol.y, 3);
    }

    #[test]
    fn test_pell_d13() {
        // x² - 13y² = 1 → fundamental: x=649, y=180
        let sol = pell_solve(13);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert_eq!(sol.x, 649);
        assert_eq!(sol.y, 180);
    }

    #[test]
    fn test_pell_perfect_square() {
        assert!(pell_solve(4).is_none());
    }

    #[test]
    fn test_pell_negative_d2() {
        // x² - 2y² = -1 → x=1, y=1
        let sol = pell_solve_negative(2);
        assert!(sol.is_some());
        let sol = sol.unwrap();
        assert_eq!(sol.x, 1);
        assert_eq!(sol.y, 1);
    }

    #[test]
    fn test_pell_nth_solution() {
        let fund = PellSolution { x: 3, y: 2 };
        let sol2 = pell_nth_solution(2, &fund, 2);
        assert_eq!(sol2.x, 17);
        assert_eq!(sol2.y, 12);
        assert_eq!(sol2.x as u128 * sol2.x as u128 - 2 * sol2.y as u128 * sol2.y as u128, 1);
    }
}
