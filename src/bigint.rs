use crate::{
    impl_arithmetic_op, impl_divlike_op,
    shared::{Bounded, DEFAULT_LIMBS},
};
use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod, Uint, Wrapping, Zero};
use std::{
    clone::Clone,
    convert::From,
    fmt,
    ops::{Add, Div, Mul, Rem, Sub},
};

/// Simple wrapper to abstract details away from crypto_bigint library.
/// We simply want to be able to:
/// do arithmetics with BigInt,
/// display it with println! for debug purposes,
/// create a BigInt from a regular integer,
/// generate a random BigInt and
/// change the size of a BigInt
#[derive(PartialEq, PartialOrd, Clone, Copy, Debug)]
pub struct BigInt<const L: usize = DEFAULT_LIMBS>(pub Wrapping<Uint<L>>);

impl<const L: usize> Bounded for BigInt<L> {
    const L: usize = L;
}

impl<const L: usize> BigInt<L> {
    /// Creates a random BigInt modulo `modulus`
    pub fn random_mod(modulus: &BigInt<L>) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(Uint::<L>::random_mod(
            &mut OsRng,
            &NonZero::new(modulus.to_uint()).unwrap(),
        )))
    }

    /// Wraps a Uint into a BigInt
    pub fn new(n: Uint<L>) -> BigInt<L> {
        BigInt(Wrapping(n))
    }

    /// Returns underlying Uint
    pub fn to_uint(self) -> Uint<L> {
        self.0 .0
    }

    /// Returns `true` if it is zero
    pub fn is_zero(&self) -> bool {
        self.0.is_zero().into()
    }

    /// Resizes a BigInt
    pub fn resize<const LNEW: usize>(&self) -> BigInt<LNEW> {
        BigInt::new(self.to_uint().resize::<LNEW>())
    }

    /// Computes square root of BigInt
    pub fn sqrt(self) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(self.to_uint().sqrt_vartime()))
    }

    /// Computes gcd of two BigInt, the good-old Euclid way
    pub fn gcd(b1: &BigInt<L>, b2: &BigInt<L>) -> BigInt<L> {
        if b1 < b2 {
            return Self::gcd(b2, b1);
        }
        let (mut x0, mut x1) = (b1.clone(), b2.clone());
        while !x1.is_zero() {
            (x1, x0) = (x0 % x1, x1);
        }
        x0
    }
}

/// Creates a BigInt from a regular integer
impl<const L: usize> From<u128> for BigInt<L> {
    fn from(k: u128) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(Uint::<L>::from(k)))
    }
}

/// Creates a BigInt from a regular &integer
impl<const L: usize> From<&u128> for BigInt<L> {
    fn from(k: &u128) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(Uint::<L>::from(*k)))
    }
}

/// Displays a BigInt
impl<const L: usize> fmt::Display for BigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "0x{}",
            self.to_uint().to_string().trim_start_matches('0')
        )
    }
}

impl_arithmetic_op!(Add, add, +);
impl_arithmetic_op!(Sub, sub, -);
impl_arithmetic_op!(Mul, mul, *);

impl_divlike_op!(Div, div, /);
impl_divlike_op!(Rem, rem, %);

#[cfg(test)]
mod tests {
    use super::*;

    const L: usize = DEFAULT_LIMBS;

    #[test]
    fn add_bigint() {
        fn simple_tester(a: u128, b: u128) {
            let big_a: BigInt<L> = BigInt::from(a);
            let big_b = BigInt::from(b);
            assert_eq!(
                (&big_a + &big_b).to_uint(),
                big_a.to_uint().wrapping_add(&big_b.to_uint())
            )
        }

        simple_tester(0, 1);
        simple_tester(129812, 92373829187);
    }
}
