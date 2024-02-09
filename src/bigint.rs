use crate::shared::DEFAULT_LIMBS;
use crypto_bigint::{
    rand_core::OsRng,
    subtle::{Choice, ConditionallySelectable, ConstantTimeEq, ConstantTimeLess},
    Checked, NonZero, RandomMod, Uint, Wrapping, Zero,
};
use std::{clone::Clone, convert::From, fmt};

/// A trait that define a big int interface. We need to do basic arithmetic operations with them,
/// computing greater common divisor, square root, generate a random int mod `modulus`, cast a u128
/// into a big int.
pub trait BigIntTrait: PartialEq + PartialOrd + Clone + fmt::Display {
    fn add<'a, 'b>(&'a self, other: &'b Self) -> Self;
    fn sub<'a, 'b>(&'a self, other: &'b Self) -> Self;
    fn mul<'a, 'b>(&'a self, other: &'b Self) -> Self;
    fn div<'a, 'b>(&'a self, other: &'b Self) -> Self;
    fn rem<'a, 'b>(&'a self, other: &'b Self) -> Self;
    fn gcd<'a, 'b>(&'a self, other: &'b Self) -> Self;
    fn sqrt<'a>(&'a self) -> Self;
    fn random_mod(modulus: &Self) -> Self;
    fn from_u128(n: u128) -> Self;
    fn is_zero<'a>(&'a self) -> Choice;
}

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct WrappingCryptoBigInt<const L: usize = DEFAULT_LIMBS>(pub Wrapping<Uint<L>>);
// #[derive(PartialEq, PartialOrd, Clone, Debug)]
#[derive(Clone, Debug)]
pub struct CheckedCryptoBigInt<const L: usize = DEFAULT_LIMBS>(pub Checked<Uint<L>>);

/// Displays a Wrapping big int
impl<const L: usize> fmt::Display for WrappingCryptoBigInt<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "0x{}", self.0 .0.to_string().trim_start_matches('0'))
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
        } else if !self_is_some && !other_is_some {
            true
        } else {
            false
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

impl<const L: usize> BigIntTrait for WrappingCryptoBigInt<L> {
    fn add<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 + other.0)
    }
    fn sub<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 - other.0)
    }
    fn mul<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 * other.0)
    }
    fn div<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 / NonZero::new(other.0 .0).unwrap())
    }
    fn rem<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 % NonZero::new(other.0 .0).unwrap())
    }
    /// Computes gcd of two &BigInt, the good-old Euclid way
    fn gcd<'a, 'b>(&'a self, other: &'b Self) -> Self {
        if self < other {
            return other.gcd(self);
        }
        let (mut x0, mut x1) = (self.clone(), other.clone());
        while !<Choice as Into<bool>>::into(x1.is_zero()) {
            (x1, x0) = (x0.rem(&x1), x1);
        }
        x0
    }
    fn sqrt<'a>(&'a self) -> Self {
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
    fn is_zero<'a>(&'a self) -> Choice {
        self.0.is_zero()
    }
}

impl<const L: usize> BigIntTrait for CheckedCryptoBigInt<L> {
    fn add<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 + other.0)
    }
    fn sub<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 - other.0)
    }
    fn mul<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(self.0 * other.0)
    }
    fn div<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(Checked(
            self.0
                 .0
                .and_then(|n| other.0 .0.and_then(|m| n.checked_div(&m))),
        ))
    }
    fn rem<'a, 'b>(&'a self, other: &'b Self) -> Self {
        Self(Checked(
            self.0
                 .0
                .and_then(|n| other.0 .0.and_then(|m| n.checked_rem(&m))),
        ))
    }
    /// Computes gcd of two &BigInt, the good-old Euclid way
    fn gcd<'a, 'b>(&'a self, other: &'b Self) -> Self {
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
    fn sqrt<'a>(&'a self) -> Self {
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
    fn is_zero<'a>(&'a self) -> Choice {
        self.0.ct_eq(&Self::from_u128(0).0)
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
}
