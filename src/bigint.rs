use crate::shared::{Bounded, DEFAULT_LIMBS};
use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod, Uint, Wrapping, Zero, U128};
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
/// - change the size of BigInt
#[derive(PartialEq, PartialOrd)]
pub struct BigInt<const L: usize = DEFAULT_LIMBS>(pub Wrapping<Uint<L>>);

impl<const L: usize> Bounded for BigInt<L> {
    const L: usize = L;
}

impl<const L: usize> BigInt<L> {
    /// create a random BigInt modulo `modulus`
    pub fn random_mod(modulus: &BigInt<L>) -> BigInt<L> {
        BigInt(Wrapping(Uint::<L>::random_mod(
            &mut OsRng,
            &NonZero::new(modulus.0 .0).unwrap(),
        )))
    }

    /// wrap a Uint into a BigInt
    pub fn new(n: Uint<L>) -> BigInt<L> {
        BigInt(Wrapping(n))
    }

    pub fn resize<const Lnew: usize>(&self) -> BigInt<Lnew> {
        BigInt::new(self.0 .0.resize::<Lnew>())
    }

    /// compute square root of BigInt
    pub fn sqrt(self) -> BigInt<L> {
        BigInt(Wrapping(self.0 .0.sqrt_vartime()))
    }

    /// compute gcd of two BigInt, the good-old Euclid way
    pub fn gcd(b1: &BigInt<L>, b2: &BigInt<L>) -> BigInt<L> {
        if b1 < b2 {
            return Self::gcd(b2, b1);
        }
        let (mut x0, mut x1) = (b1.clone(), b2.clone());
        while x1.0 .0.is_zero().unwrap_u8() == 0 {
            (x1, x0) = (&x0 % &x1, x1);
        }
        return x0;
    }
}
// clone a BgInt
impl<const L: usize> Clone for BigInt<L> {
    fn clone(&self) -> Self {
        BigInt::new(self.0 .0.clone())
    }
}

// creates a BigInt from a regular integer
impl<const L: usize> From<u128> for BigInt<L> {
    fn from(k: u128) -> BigInt<L> {
        BigInt(Wrapping(Uint::<L>::from(k)))
    }
}

impl<const L: usize> From<&u128> for BigInt<L> {
    fn from(k: &u128) -> BigInt<L> {
        BigInt(Wrapping(Uint::<L>::from(k.clone())))
    }
}

// display a BigInt
impl<const L: usize> fmt::Display for BigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0 .0.to_string())
    }
}

// implement add, sub, mul, div, rem for &BigInt
impl<'a, 'b, const L: usize> Add<&'b BigInt<L>> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn add(self, other: &'b BigInt<L>) -> BigInt<L> {
        let (b1, b2) = (self.0, other.0);
        BigInt(b1 + b2)
    }
}
impl<'a, 'b, const L: usize> Sub<&'b BigInt<L>> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn sub(self, other: &'b BigInt<L>) -> BigInt<L> {
        let (b1, b2) = (self.0, other.0);
        BigInt(b1 - b2)
    }
}

/// macro to implement arithmetic operations for many sizes
/// `$op` is the arithmetic operation
/// `$struct` is the struct that gets implemented
/// `$l1, $l2, $l3` are some `const` that must be known at compile time
macro_rules! impl_arithmetics {
    (+, $struct: ident, $l1:expr, $l2: expr, $l3: expr) => {
        impl_arithmetics!(Add, add, +, $struct, $l1, $l2, $l3);
    };

    (-, $struct: ident, $l1:expr, $l2: expr, $l3: expr) => {
        impl_arithmetics!(Sub, sub, -, $struct, $l1, $l2, $l3);
    };

    (*, $struct: ident, $l1:expr, $l2: expr, $l3: expr) => {
        impl_arithmetics!(Mul, mul, *, $struct, $l1, $l2, $l3);
    };

    (/, $struct: ident, $l1:expr, $l2: expr, $l3: expr) => {
        impl_arithmetics!(Div, div, /, $struct, $l1, $l2, $l3);
    };

    ($trait: ident, $function: ident, $op: tt, $struct: ident, $l1:expr, $l2: expr, $l3: expr) => {
        impl<'a, 'b> $trait<&'b $struct<$l2>> for &'a $struct<$l1>
        where
            $struct<$l1>: Sized,
            $struct<$l2>: Sized,
            $struct<$l3>: Sized,
        {
            type Output = $struct<$l3>;
            fn $function(self, other: &'b $struct<$l2>) -> $struct<$l3> {
                $struct::<$l3>(self.0 $op other.0)
            }
        }
        impl $trait<$struct<$l2>> for $struct<$l1>
        where
            $struct<$l1>: Sized,
            $struct<$l2>: Sized,
            $struct<$l3>: Sized,
        {
            type Output = $struct<$l3>;
            fn $function(self, other: $struct<$l2>) -> $struct<$l3> {
                &self $op &other
            }
        }
    };
}

// TODO: repeat with all possible constants
impl_arithmetics!(+, BigInt, DEFAULT_LIMBS, DEFAULT_LIMBS , DEFAULT_LIMBS);
impl_arithmetics!(-, BigInt, DEFAULT_LIMBS, DEFAULT_LIMBS , DEFAULT_LIMBS);
impl_arithmetics!(/, BigInt, DEFAULT_LIMBS, DEFAULT_LIMBS , DEFAULT_LIMBS);
impl_arithmetics!(*, BigInt, DEFAULT_LIMBS, DEFAULT_LIMBS , {2*DEFAULT_LIMBS});
impl_arithmetics!(*, BigInt, DEFAULT_LIMBS, { 2 * DEFAULT_LIMBS }, {
    3 * DEFAULT_LIMBS
});

impl<'a, 'b, const L: usize> Div<&'b BigInt<L>> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn div(self, other: &'b BigInt<L>) -> BigInt<L> {
        let (b1, b2) = (self.0, NonZero::new(other.0 .0).unwrap());
        BigInt(b1 / b2)
    }
}

impl<'a, 'b, const L: usize> Rem<&'b BigInt<L>> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn rem(self, other: &'b BigInt<L>) -> BigInt<L> {
        let (b1, b2) = (self.0, NonZero::new(other.0 .0).unwrap());
        BigInt(b1 % b2)
    }
}

// implement operations for BigInt
impl<const L: usize> Add<BigInt<L>> for BigInt<L> {
    type Output = BigInt<L>;
    fn add(self, other: BigInt<L>) -> BigInt<L> {
        &self + &other
    }
}

impl<const L: usize> Sub<BigInt<L>> for BigInt<L> {
    type Output = BigInt<L>;
    fn sub(self, other: BigInt<L>) -> BigInt<L> {
        &self - &other
    }
}

impl<const L1: usize, const L2: usize> Mul<BigInt<L2>> for BigInt<L1>
where
    BigInt<{ L1 + L2 }>: Sized,
{
    type Output = BigInt<{ L1 + L2 }>;
    fn mul(self, other: BigInt<L2>) -> BigInt<{ L1 + L2 }> {
        &self * &other
    }
}

impl<const L: usize> Div<BigInt<L>> for BigInt<L> {
    type Output = BigInt<L>;
    fn div(self, other: BigInt<L>) -> BigInt<L> {
        &self / &other
    }
}

impl<const L: usize> Rem<BigInt<L>> for BigInt<L> {
    type Output = BigInt<L>;
    fn rem(self, other: BigInt<L>) -> BigInt<L> {
        &self % &other
    }
}

// implement arithmetic operations with u128
impl<'a, const L: usize> Add<u128> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn add(self, other: u128) -> BigInt<L> {
        self + &BigInt::from(other)
    }
}
impl<'a, const L: usize> Sub<u128> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn sub(self, other: u128) -> BigInt<L> {
        self - &BigInt::from(other)
    }
}

impl<'a, const L: usize> Mul<u128> for &'a BigInt<L>
where
    BigInt<{ L + U128::LIMBS }>: Sized,
    BigInt<{ U128::LIMBS }>: Sized,
{
    type Output = BigInt<{ L + U128::LIMBS }>;
    fn mul(self, other: u128) -> BigInt<{ L + U128::LIMBS }> {
        self * &BigInt::<{ U128::LIMBS }>::from(other)
    }
}
impl<'a, const L: usize> Div<u128> for &'a BigInt<L> {
    type Output = BigInt<L>;
    fn div(self, other: u128) -> BigInt<L> {
        self / &BigInt::from(other)
    }
}
