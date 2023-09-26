#[macro_use]
extern crate impl_ops;

use crate::shared::{Bounded, DEFAULT_LIMBS};
use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod, Uint, Wrapping, Zero, U128};
use std::{clone::Clone, convert::From, fmt, ops};

/// Simple wrapper to abstract details away from crypto_bigint library.
/// We simply want to be able to:
/// - do arithmetics with BigInt
/// - display it with println! for debug purposes
/// - create a BigInt from a regular integer
/// - generate a random BigInt
/// - change the size of BigInt
#[derive(PartialEq, PartialOrd)]
pub struct BigInt<const L: usize = DEFAULT_LIMBS>(pub Uint<L>);

impl<const L: usize> Bounded for BigInt<L> {
    const L: usize = L;
}

impl<const L: usize> BigInt<L> {
    /// create a random BigInt modulo `modulus`
    pub fn random_mod(modulus: &BigInt<L>) -> BigInt<L> {
        BigInt(Uint::<L>::random_mod(
            &mut OsRng,
            &NonZero::new(modulus.0).unwrap(),
        ))
    }

    /// wrap a Uint into a BigInt
    pub fn new(n: Uint<L>) -> BigInt<L> {
        BigInt(n)
    }

    /// resize a BigInt
    pub fn resize<const Lnew: usize>(&self) -> BigInt<Lnew> {
        BigInt::new(self.0.resize::<Lnew>())
    }

    /// compute square root of BigInt
    pub fn sqrt(self) -> BigInt<L> {
        BigInt(self.0.sqrt_vartime())
    }

    /// compute gcd of two BigInt, the good-old Euclid way
    pub fn gcd(b1: &BigInt<L>, b2: &BigInt<L>) -> BigInt<L> {
        if b1 < b2 {
            return Self::gcd(b2, b1);
        }
        let (mut x0, mut x1) = (b1.clone(), b2.clone());
        while x1.0.is_zero().unwrap_u8() == 0 {
            (x1, x0) = (&x0 % &x1, x1);
        }
        return x0;
    }
}
/// clone a BigInt
impl<const L: usize> Clone for BigInt<L> {
    fn clone(&self) -> Self {
        BigInt::new(self.0.clone())
    }
}

/// creates a BigInt from a regular integer
impl From<u128> for BigInt<{ U128::LIMBS }> {
    fn from(k: u128) -> BigInt<{ U128::LIMBS }> {
        BigInt(Uint::<{ U128::LIMBS }>::from(k))
    }
}

/// creates a BigInt from a regular &integer
impl<const L: usize> From<&u128> for BigInt<L> {
    fn from(k: &u128) -> BigInt<L> {
        BigInt(Uint::<L>::from(k.clone()))
    }
}

/// display a BigInt
impl<const L: usize> fmt::Display for BigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0.to_string())
    }
}

/// macro to implement arithmetic operations for many sizes
/// `$struct` is the struct that gets implemented
/// `$l1, $l2, $l3` are some `const` that must be known at compile time
macro_rules! impl_sized_ops {
    (+, $type1: path, $type2: path, $type3:path, $size: expr) => {
        impl_op_ex!( + | a: $type1, b: $type2 | -> $type3
            {BigInt((Wrapping(a.0.resize::<$size>()) + Wrapping(b.0.resize::<$size>())).0) }
        );
    };

    (-, $type1: path, $type2: path, $type3:path, $size: expr) => {
        impl_op_ex!( - | a: $type1, b: $type2 | -> $type3
            { BigInt((Wrapping(a.0.resize::<$size>()) - Wrapping(b.0.resize::<$size>())).0) }
        );
    };

    (/, $type1: path, $type2: path, $type3:path, $size: expr) => {
        impl_op_ex!( / | a: $type1, b: $type2 | -> $type3
            {BigInt( a.0 / NonZero::new(b.0).unwrap() )}
        );
    };

    (%, $type1: path, $type2: path, $type3:path, $size: expr) => {
        impl_op_ex!( % | a: $type1, b: $type2 | -> $type3
            {BigInt( a.0 % NonZero::new(b.0).unwrap() )}
        );
    };

    (*, $type1: path, $type2: path, $type3:path, $size: expr) => {
        impl_op_ex!( * | a: $type1, b: $type2 | -> $type3
            {BigInt::<$size>( a.0 * b.0 )}
        );
    };
}

macro_rules! impl_bigint_sized_ops {
    ($mult1: expr, $mult2: expr, $mult3: expr) => {
        impl_sized_ops!(
            *,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult2 * U128::LIMBS }>,
            BigInt<{ $mult3 * U128::LIMBS }>,
            { $mult3 * U128::LIMBS }
        );
    };
    ($mult1: expr) => {
        impl_sized_ops!(
            +,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            { $mult1 * U128::LIMBS }
        );
        impl_sized_ops!(
            -,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            { $mult1 * U128::LIMBS }
        );
        impl_sized_ops!(
            /,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            { $mult1 * U128::LIMBS }
        );
        impl_sized_ops!(
            %,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            BigInt<{ $mult1 * U128::LIMBS }>,
            { $mult1 * U128::LIMBS }
        );
    };
}

// TODO: repeat with all possible constants
impl_bigint_sized_ops!(1);
impl_bigint_sized_ops!(2);
impl_bigint_sized_ops!(3);
impl_bigint_sized_ops!(4);
impl_bigint_sized_ops!(5);
impl_bigint_sized_ops!(6);
impl_bigint_sized_ops!(7);
impl_bigint_sized_ops!(8);

impl_bigint_sized_ops!(1, 1, 2);
impl_bigint_sized_ops!(1, 2, 3);
impl_bigint_sized_ops!(1, 3, 4);
impl_bigint_sized_ops!(1, 4, 5);
impl_bigint_sized_ops!(1, 5, 6);
impl_bigint_sized_ops!(1, 6, 7);
impl_bigint_sized_ops!(1, 7, 8);
impl_bigint_sized_ops!(1, 8, 9);
impl_bigint_sized_ops!(2, 1, 3);
impl_bigint_sized_ops!(2, 2, 4);
impl_bigint_sized_ops!(2, 3, 5);
impl_bigint_sized_ops!(2, 4, 6);
impl_bigint_sized_ops!(2, 5, 7);
impl_bigint_sized_ops!(2, 6, 8);
impl_bigint_sized_ops!(2, 7, 9);
impl_bigint_sized_ops!(2, 8, 10);
impl_bigint_sized_ops!(3, 1, 4);
impl_bigint_sized_ops!(3, 2, 5);
impl_bigint_sized_ops!(3, 3, 6);
impl_bigint_sized_ops!(3, 4, 7);
impl_bigint_sized_ops!(3, 5, 8);
impl_bigint_sized_ops!(3, 6, 9);
impl_bigint_sized_ops!(3, 7, 10);
impl_bigint_sized_ops!(3, 8, 11);
impl_bigint_sized_ops!(4, 1, 5);
impl_bigint_sized_ops!(4, 2, 6);
impl_bigint_sized_ops!(4, 3, 7);
impl_bigint_sized_ops!(4, 4, 8);
impl_bigint_sized_ops!(4, 5, 9);
impl_bigint_sized_ops!(4, 6, 10);
impl_bigint_sized_ops!(4, 7, 11);
impl_bigint_sized_ops!(4, 8, 12);
impl_bigint_sized_ops!(5, 1, 6);
impl_bigint_sized_ops!(5, 2, 7);
impl_bigint_sized_ops!(5, 3, 8);
impl_bigint_sized_ops!(5, 4, 9);
impl_bigint_sized_ops!(5, 5, 10);
impl_bigint_sized_ops!(5, 6, 11);
impl_bigint_sized_ops!(5, 7, 12);
impl_bigint_sized_ops!(5, 8, 13);
impl_bigint_sized_ops!(6, 1, 7);
impl_bigint_sized_ops!(6, 2, 8);
impl_bigint_sized_ops!(6, 3, 9);
impl_bigint_sized_ops!(6, 4, 10);
impl_bigint_sized_ops!(6, 5, 11);
impl_bigint_sized_ops!(6, 6, 12);
impl_bigint_sized_ops!(6, 7, 13);
impl_bigint_sized_ops!(6, 8, 14);
impl_bigint_sized_ops!(7, 1, 8);
impl_bigint_sized_ops!(7, 2, 9);
impl_bigint_sized_ops!(7, 3, 10);
impl_bigint_sized_ops!(7, 4, 11);
impl_bigint_sized_ops!(7, 5, 12);
impl_bigint_sized_ops!(7, 6, 13);
impl_bigint_sized_ops!(7, 7, 14);
impl_bigint_sized_ops!(7, 8, 15);
impl_bigint_sized_ops!(8, 1, 9);
impl_bigint_sized_ops!(8, 2, 10);
impl_bigint_sized_ops!(8, 3, 11);
impl_bigint_sized_ops!(8, 4, 12);
impl_bigint_sized_ops!(8, 5, 13);
impl_bigint_sized_ops!(8, 6, 14);
impl_bigint_sized_ops!(8, 7, 15);
impl_bigint_sized_ops!(8, 8, 16);

#[cfg(test)]
mod tests {
    use super::*;

    const L: usize = DEFAULT_LIMBS;

    #[test]
    fn add_bigint() {
        fn simple_tester(a: u128, b: u128) {
            let big_a = BigInt::from(a);
            let big_b = BigInt::from(b);
            assert_eq!((big_a + big_b).0, big_a.0.wrapping_add(big_b.0))
        }

        simple_tester(0, 1);
        simple_tester(129812, 92373829187);
    }
}
