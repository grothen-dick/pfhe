use crate::shared::{Bounded, DEFAULT_LIMBS};
use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod, Uint, Wrapping, Zero};
use std::{
    clone::Clone,
    convert::From,
    fmt,
    ops::{Add, Div, Mul, Rem, Sub},
};

/// Simple wrapper to abstract details away from crypto_bigint library.
/// We simply want to be able to:
/// - do arithmetics with BigInt
/// - display it with println! for debug purposes
/// - create a BigInt from a regular integer
/// - generate a random BigInt
/// - change the size of a BigInt
#[derive(PartialEq, PartialOrd, Clone, Copy)]
pub struct BigInt<const L: usize = DEFAULT_LIMBS>(pub Wrapping<Uint<L>>);

impl<const L: usize> Bounded for BigInt<L> {
    const L: usize = L;
}

impl<const L: usize> BigInt<L> {
    /// create a random BigInt modulo `modulus`
    pub fn random_mod(modulus: &BigInt<L>) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(Uint::<L>::random_mod(
            &mut OsRng,
            &NonZero::new(modulus.0 .0).unwrap(),
        )))
    }

    /// wrap a Uint into a BigInt
    pub fn new(n: Uint<L>) -> BigInt<L> {
        BigInt(Wrapping(n))
    }

    /// resize a BigInt
    pub fn resize<const LNEW: usize>(&self) -> BigInt<LNEW> {
        BigInt::new(self.0 .0.resize::<LNEW>())
    }

    /// compute square root of BigInt
    pub fn sqrt(self) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(self.0 .0.sqrt_vartime()))
    }

    /// compute gcd of two BigInt, the good-old Euclid way
    pub fn gcd(b1: &BigInt<L>, b2: &BigInt<L>) -> BigInt<L> {
        if b1 < b2 {
            return Self::gcd(b2, b1);
        }
        let (mut x0, mut x1) = (*b1, *b2);
        while x1.0.is_zero().unwrap_u8() == 0 {
            (x1, x0) = (x0 % x1, x1);
        }
        x0
    }
}

/// creates a BigInt from a regular integer
impl<const L: usize> From<u128> for BigInt<L> {
    fn from(k: u128) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(Uint::<L>::from(k)))
    }
}

/// creates a BigInt from a regular &integer
impl<const L: usize> From<&u128> for BigInt<L> {
    fn from(k: &u128) -> BigInt<L> {
        BigInt(Wrapping::<Uint<L>>(Uint::<L>::from(*k)))
    }
}

/// display a BigInt
impl<const L: usize> fmt::Display for BigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

macro_rules! impl_arithmetic_op {
    ($trait: ident, $function: ident, $op: tt) => {
        impl<'a, 'b, const L: usize> $trait<&'b BigInt<L>> for &'a BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: &'b BigInt<L>) -> BigInt<L> {
                BigInt(self.0 $op other.0)
            }
        }

        impl<const L: usize> $trait<BigInt<L>> for BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: BigInt<L>) -> BigInt<L> {
                &self $op &other
            }
        }

    };
}

impl_arithmetic_op!(Add, add, +);
impl_arithmetic_op!(Sub, sub, -);
impl_arithmetic_op!(Mul, mul, *);

macro_rules! impl_divlike_op {
    ($trait: ident, $function: ident, $op: tt) => {
        impl<'a, 'b, const L: usize> $trait<&'b BigInt<L>> for &'a BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: &'b BigInt<L>) -> BigInt<L> {
                BigInt(self.0 $op NonZero::new(other.0.0).unwrap())
            }
        }

        impl<const L: usize> $trait<BigInt<L>> for BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: BigInt<L>) -> BigInt<L> {
                &self $op &other
            }
        }

    };
}

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
            assert_eq!((&big_a + &big_b).0 .0, big_a.0 .0.wrapping_add(&big_b.0 .0))
        }

        simple_tester(0, 1);
        simple_tester(129812, 92373829187);
    }
}
