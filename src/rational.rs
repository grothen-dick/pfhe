use super::{
    fmt,
    hensel_code::HenselCode,
    ops::{Add, Mul},
    Clone,
};
use crypto_bigint::subtle::Choice;

use crate::bigint::BigIntTrait;

#[derive(Clone, PartialEq, Debug)]
pub struct Rational<T: BigIntTrait> {
    pub num: T,
    pub denom: T,
}

/// Basic functionalities: simplify common factors, resize
impl<T: BigIntTrait> Rational<T> {
    pub fn reduce(&self) -> Self {
        let gcd = &self.num.gcd(&self.denom);
        let num = self.num.div(gcd);
        let denom = self.denom.div(gcd);
        Rational::<T> { num, denom }
    }
}

/// Pretty-prints Rational
impl<T: BigIntTrait> fmt::Display for Rational<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// Adds two &Rational
impl<'a, 'b, T: BigIntTrait> Add<&'b Rational<T>> for &'a Rational<T> {
    type Output = Rational<T>;
    fn add(self, other: &'b Rational<T>) -> Rational<T> {
        let num1 = self.num.mul(&other.denom);
        let num2 = self.denom.mul(&other.num);
        let num = num1.add(&num2);
        let denom = self.denom.mul(&other.denom);
        Rational::<T> { num, denom }.reduce()
    }
}

/// Adds two Rational
impl<T: BigIntTrait> Add<Rational<T>> for Rational<T> {
    type Output = Rational<T>;
    fn add(self, other: Rational<T>) -> Rational<T> {
        &self + &other
    }
}

/// Multiplies two &Rational
impl<'a, 'b, T: BigIntTrait> Mul<&'b Rational<T>> for &'a Rational<T> {
    type Output = Rational<T>;
    fn mul(self, other: &'b Rational<T>) -> Rational<T> {
        Rational::<T> {
            num: self.num.mul(&other.num),
            denom: self.denom.mul(&other.denom),
        }
        .reduce()
    }
}

/// Multiplies two Rational
impl<T: BigIntTrait> Mul<Rational<T>> for Rational<T> {
    type Output = Rational<T>;
    fn mul(self, other: Rational<T>) -> Rational<T> {
        &self * &other
    }
}

/// Given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)), returns a rational
/// num/denom where:
///  i)   0 <= num   <= n_max,
///  ii)  0 <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
impl<T: BigIntTrait> From<&HenselCode<T>> for Rational<T> {
    fn from(hc: &HenselCode<T>) -> Self {
        let n_max = hc
            .modulus
            .sub(&T::from_u128(1))
            .div(&T::from_u128(2))
            .sqrt();

        let (mut x0, mut x1) = (hc.modulus.clone(), hc.res.clone());

        if x1.is_zero().into() {
            return Rational::<T> {
                num: T::from_u128(0),
                denom: T::from_u128(1),
            };
        }
        // perform (modified) extended euclidean algorithm on (g, n % g)
        let (mut y0, mut y1) = (T::from_u128(0), T::from_u128(1));
        while (x0 > n_max) && !<Choice as Into<bool>>::into(x1.is_zero()) {
            let q = x0.div(&x1);
            let (new_x0, new_x1) = (x1.clone(), x0.sub(&q.mul(&x1)));
            let (new_y0, new_y1) = (y1.clone(), y0.sub(&q.mul(&y1)));
            (x0, x1) = (new_x0, new_x1);
            (y0, y1) = (new_y0, new_y1);
        }

        Rational::<T> { num: x0, denom: y0 }
    }
}

#[cfg(test)]
mod tests {
    use super::Rational;
    use crate::bigint::{BigIntTrait, WrappingCryptoBigInt};
    // use crate::shared::DEFAULT_LIMBS;

    // const L: usize = DEFAULT_LIMBS;
    type T = WrappingCryptoBigInt;

    #[test]
    fn adds_rationals() {
        fn simple_tester(r1: &Rational<T>, r2: &Rational<T>) -> () {
            let sum = r1 + r2;
            assert_eq!(sum.denom, (r1.denom.mul(&r2.denom)));
            assert_eq!(sum.num, r1.denom.mul(&r2.num).add(&r2.denom.mul(&r1.num)));
            println!("{} + {} = {}", r1, r2, sum);
        }
        let from_u128 = <T as BigIntTrait>::from_u128;

        // integer addition
        let r11 = Rational::<T> {
            num: from_u128(3),
            denom: from_u128(1),
        };
        let r21 = Rational::<T> {
            num: from_u128(1312),
            denom: from_u128(1),
        };
        simple_tester(&r11, &r21);

        // rational addition
        let r12 = Rational::<T> {
            num: from_u128(1337),
            denom: from_u128(41),
        };
        let r22 = Rational::<T> {
            num: from_u128(57),
            denom: from_u128(982),
        };
        simple_tester(&r12, &r22);
    }
}
