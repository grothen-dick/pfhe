extern crate num_bigint_dig;
extern crate rand;

use num_bigint_dig::{BigInt, RandPrime, RandomBits};
use rand::{thread_rng, Rng};

use crate::shared::DEFAULT_LIMBS;

use crypto_bigint::{
    rand_core::OsRng,
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, ConstantTimeLess},
    Checked, NonZero, RandomMod, Uint, Wrapping, Zero,
};

use crypto_primes::generate_prime as crypto_primes_generate;

use std::{clone::Clone, fmt};

/// A trait that define a big int interface. We need to do basic arithmetic operations with them,
/// computing greater common divisor, square root, generate a random int mod `modulus`, cast a u128
/// into a big int.
pub trait BigIntTrait: PartialEq + PartialOrd + Clone + fmt::Display + fmt::Debug {
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn pow(&self, other: u128) -> Self {
        let mut exponent = other;
        if exponent == 0 {
            return Self::one();
        }
        let mut x = self.clone();
        let mut y = Self::one();
        while exponent > 1 {
            if exponent % 2 == 1 {
                y = x.mul(&y);
                exponent -= 1;
            }
            x = x.mul(&x);
            exponent /= 2;
        }
        x.mul(&y)
    }
    fn div(&self, other: &Self) -> Self;
    fn rem(&self, other: &Self) -> Self;
    fn gcd(&self, other: &Self) -> Self;
    fn sqrt(&self) -> Self;
    fn random_mod(modulus: &Self) -> Self;
    fn from_u128(n: u128) -> Self;
    fn is_zero(&self) -> Choice;
    fn generate_prime(bit_length: Option<usize>) -> Self;
    fn zero() -> Self {
        Self::from_u128(0)
    }
    fn one() -> Self {
        Self::from_u128(1)
    }
}

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct WrappingCryptoBigInt<const L: usize = DEFAULT_LIMBS>(pub Wrapping<Uint<L>>);
// #[derive(PartialEq, PartialOrd, Clone, Debug)]
#[derive(Clone, Debug)]
pub struct CheckedCryptoBigInt<const L: usize = DEFAULT_LIMBS>(pub Checked<Uint<L>>);

/// Displays a Wrapping big int
impl<const L: usize> fmt::Display for WrappingCryptoBigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.0 .0.eq(&Uint::<L>::ZERO) {
            write!(f, "0x0")
        } else {
            write!(f, "0x{}", self.0 .0.to_string().trim_start_matches('0'))
        }
    }
}

/// Displays a Checked big int
impl<const L: usize> fmt::Display for CheckedCryptoBigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let b: bool = self.0 .0.is_some().into();
        if b {
            let n = self.0 .0.unwrap();
            write!(f, "{:?}", n.to_string().trim_start_matches('0'))
        } else {
            write!(f, "(NONE VALUE)")
        }
    }
}

/// Checks wether two Checked big int are equal
impl<const L: usize> std::cmp::PartialEq for CheckedCryptoBigInt<L> {
    fn eq(&self, other: &Self) -> bool {
        let (self_is_some, other_is_some): (bool, bool) =
            (self.0 .0.is_some().into(), other.0 .0.is_some().into());
        if self_is_some && other_is_some {
            self.0 .0.unwrap() == other.0 .0.unwrap()
        } else {
            !self_is_some && !other_is_some
        }
    }
}

/// Compares two Checked big int
impl<const L: usize> std::cmp::PartialOrd for CheckedCryptoBigInt<L> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let (self_is_some, other_is_some): (bool, bool) =
            (self.0 .0.is_some().into(), other.0 .0.is_some().into());
        if self_is_some && other_is_some {
            self.0 .0.unwrap().partial_cmp(&other.0 .0.unwrap())
        } else {
            None
        }
    }
}

impl BigIntTrait for BigInt {
    fn add(&self, other: &Self) -> Self {
        self + other
    }
    fn sub(&self, other: &Self) -> Self {
        self - other
    }
    fn mul(&self, other: &Self) -> Self {
        self * other
    }
    fn div(&self, other: &Self) -> Self {
        self / other
    }
    fn rem(&self, other: &Self) -> Self {
        let result = self % other;
        if result < Self::from(0) {
            result + other
        } else {
            result
        }
    }

    /// Computes gcd of two &BigInt, the good-old Euclid way
    fn gcd(&self, other: &Self) -> Self {
        if self < other {
            return other.gcd(self);
        }
        let (mut x0, mut x1) = (self.clone(), other.clone());
        while !<Choice as Into<bool>>::into(x1.is_zero()) {
            (x1, x0) = (x0.rem(&x1), x1);
        }
        x0
    }

    fn sqrt(&self) -> Self {
        self.sqrt()
    }

    fn from_u128(n: u128) -> Self {
        Self::from(n)
    }

    // Hacky implementation because of the return type defined in crypto_bigint
    fn is_zero(&self) -> Choice {
        if *self == Self::from(0_u128) {
            Uint::<2>::from(0_u128).is_zero()
        } else {
            Uint::<2>::from(1_u128).is_zero()
        }
    }

    fn generate_prime(bit_length: Option<usize>) -> Self {
        let mut rng = thread_rng();
        if let Some(bits) = bit_length {
            Self::from(rng.gen_prime(bits))
        } else {
            // No good choice of default number of bits with num_bigint
            Self::from(rng.gen_prime(32))
        }
    }

    fn random_mod(modulus: &Self) -> Self {
        if *modulus < BigInt::from(0_u128) {
            panic!("Try to generate a random BigInt modulo a negative number")
        }
        let mut rng = thread_rng();
        let large_random_number: Self = rng.sample(RandomBits::new(modulus.bits() + 1));
        large_random_number.rem(modulus)
    }
}

impl<const L: usize> BigIntTrait for WrappingCryptoBigInt<L> {
    fn add(&self, other: &Self) -> Self {
        Self(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        Self(self.0 - other.0)
    }
    fn mul(&self, other: &Self) -> Self {
        Self(self.0 * other.0)
    }
    fn div(&self, other: &Self) -> Self {
        Self(self.0 / NonZero::new(other.0 .0).unwrap())
    }
    fn rem(&self, other: &Self) -> Self {
        Self(self.0 % NonZero::new(other.0 .0).unwrap())
    }
    /// Computes gcd of two &BigInt, the good-old Euclid way
    fn gcd(&self, other: &Self) -> Self {
        if self < other {
            return other.gcd(self);
        }
        let (mut x0, mut x1) = (self.clone(), other.clone());
        while !<Choice as Into<bool>>::into(x1.is_zero()) {
            (x1, x0) = (x0.rem(&x1), x1);
        }
        x0
    }
    fn sqrt(&self) -> Self {
        Self(Wrapping::<Uint<L>>(self.0 .0.sqrt_vartime()))
    }
    fn random_mod(modulus: &Self) -> Self {
        Self(Wrapping::<Uint<L>>(Uint::<L>::random_mod(
            &mut OsRng,
            &NonZero::new(modulus.0 .0).unwrap(),
        )))
    }
    fn from_u128(n: u128) -> Self {
        Self(Wrapping::<Uint<L>>(Uint::<L>::from(n)))
    }
    fn is_zero(&self) -> Choice {
        self.0.is_zero()
    }
    fn generate_prime(bit_length: Option<usize>) -> Self {
        WrappingCryptoBigInt(Wrapping(crypto_primes_generate::<L>(bit_length)))
    }
}

impl<const L: usize> BigIntTrait for CheckedCryptoBigInt<L> {
    fn add(&self, other: &Self) -> Self {
        Self(self.0 + other.0)
    }
    fn sub(&self, other: &Self) -> Self {
        Self(self.0 - other.0)
    }
    fn mul(&self, other: &Self) -> Self {
        Self(self.0 * other.0)
    }
    fn div(&self, other: &Self) -> Self {
        Self(Checked(
            self.0
                 .0
                .and_then(|n| other.0 .0.and_then(|m| n.checked_div(&m))),
        ))
    }
    fn rem(&self, other: &Self) -> Self {
        Self(Checked(
            self.0
                 .0
                .and_then(|n| other.0 .0.and_then(|m| n.checked_rem(&m))),
        ))
    }
    /// Computes gcd of two &BigInt, the good-old Euclid way
    fn gcd(&self, other: &Self) -> Self {
        Self(Checked(self.0 .0.and_then(|n| {
            other.0 .0.map(|m| {
                let is_n_lt_m = n.ct_lt(&m);
                let mut x0 = Uint::<L>::conditional_select(&n, &m, is_n_lt_m);
                let mut x1 = Uint::<L>::conditional_select(&m, &n, is_n_lt_m);
                while !<Choice as Into<bool>>::into(x1.is_zero()) {
                    (x1, x0) = (x0.checked_rem(&x1).unwrap(), x1);
                }
                x0
            })
        })))
    }
    fn sqrt(&self) -> Self {
        Self(Checked::<Uint<L>>(self.0 .0.map(|x| x.sqrt_vartime())))
    }
    fn random_mod(modulus: &Self) -> Self {
        Self(Checked(modulus.0 .0.map(|n| {
            Uint::<L>::random_mod(&mut OsRng, &NonZero::new(n).unwrap())
        })))
    }
    fn from_u128(n: u128) -> Self {
        Self(Checked::new(Uint::<L>::from(n)))
    }
    fn is_zero(&self) -> Choice {
        self.0.ct_eq(&Self::zero().0)
    }
    fn generate_prime(bit_length: Option<usize>) -> Self {
        CheckedCryptoBigInt(Checked::new(crypto_primes_generate::<L>(bit_length)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const L: usize = DEFAULT_LIMBS;
    type T = WrappingCryptoBigInt<L>;

    #[test]
    fn add_bigint() {
        fn simple_tester(a: u128, b: u128) {
            let big_a = T::from_u128(a);
            let big_b = T::from_u128(b);
            assert_eq!(
                (big_a.add(&big_b)).0 .0,
                big_a.0 .0.wrapping_add(&big_b.0 .0)
            )
        }

        simple_tester(0, 1);
        simple_tester(129812, 92373829187);
    }

    #[test]
    fn test_pow() {
        fn simple_tester(a: u128, b: u128) {
            let big_a = T::from_u128(a);
            assert_eq!(big_a.pow(b), T::from_u128(a.pow(b as u32)));
        }
        simple_tester(12, 4);
        simple_tester(3, 0);
        simple_tester(12, 6);
        simple_tester(12, 9);
    }
}
