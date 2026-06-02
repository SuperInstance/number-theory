//! # numtheory
//!
//! Number theory in Rust: primes, modular arithmetic,
//! Diophantine equations, continued fractions, and more.

pub mod primes;
pub mod modular;
pub mod arithmetic;
pub mod continued_fraction;
pub mod quadratic;
pub mod diophantine;
pub mod dirichlet;

pub use primes::*;
pub use modular::*;
pub use arithmetic::*;
pub use continued_fraction::*;
pub use quadratic::*;
pub use diophantine::*;
pub use dirichlet::*;
