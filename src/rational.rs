use super::{
    fmt,
    hensel_code::HenselCode,
    ops::{Add, Mul},
    Clone,
};

use crate::bigint::BigInt;
use crate::shared::DEFAULT_LIMBS;

pub struct Rational<const L: usize> {
    pub num: BigInt<L>,
    pub denom: BigInt<L>,
}

impl<const L: usize> Clone for Rational<L> {
    fn clone(&self) -> Self {
        Rational::<L> {
            num: self.num.clone(),
            denom: self.denom.clone(),
        }
    }
}

/// pretty-print Rational
impl<const L: usize> fmt::Display for Rational<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// we want to simplify a rational
impl<const L: usize> Rational<L> {
    pub fn reduce(&self) -> Self {
        let gcd = BigInt::gcd(&self.num, &self.denom);
        let num = &self.num / &gcd;
        let denom = &self.denom / &gcd;
        Rational::<L> { num, denom }
    }
}

/// macro to implement rational arithmetics for many sizes
/// `$l1, $l2, $l3` are some `const` that must be known at compile time
macro_rules! impl_rational_arithmetics {
    ($l1:expr, $l2: expr, $l3: expr) => {
        impl<'a, 'b> Add<&'b Rational<$l2>> for &'a Rational<$l1>
        where
            Rational<$l1>: Sized,
            Rational<$l2>: Sized,
            Rational<$l3>: Sized,
        {
            type Output = Rational<$l3>;
            fn add(self, other: &'b Rational<$l2>) -> Rational<$l3> {
                let (r1, r2) = (self, other);
                let num1 = &r1.num * &r2.denom;
                let num2 = &r1.denom * &r2.num;
                let num: BigInt<$l3> = num1 + num2;
                let denom: BigInt<$l3> = &r1.denom * &r2.denom;
                Rational::<$l3> { num, denom }.reduce()
            }
        }
        impl Add<Rational<$l2>> for Rational<$l1>
        where
            Rational<$l1>: Sized,
            Rational<$l2>: Sized,
            Rational<$l3>: Sized,
        {
            type Output = Rational<$l3>;
            fn add(self, other: Rational<$l2>) -> Rational<$l3> {
                &self + &other
            }
        }
        impl<'a, 'b> Mul<&'b Rational<$l2>> for &'a Rational<$l1>
        where
            Rational<$l1>: Sized,
            Rational<$l2>: Sized,
            Rational<$l3>: Sized,
        {
            type Output = Rational<$l3>;
            fn mul(self, other: &'b Rational<$l2>) -> Rational<$l3> {
                let (r1, r2) = (self, other);
                Rational::<$l3> {
                    num: &r1.num * &r2.num,
                    denom: &r1.denom * &r2.denom,
                }
                .reduce()
            }
        }
        impl Mul<Rational<$l2>> for Rational<$l1>
        where
            Rational<$l1>: Sized,
            Rational<$l2>: Sized,
            Rational<$l3>: Sized,
        {
            type Output = Rational<$l3>;
            fn mul(self, other: Rational<$l2>) -> Rational<$l3> {
                &self * &other
            }
        }
    };
}

// TODO: repeat with all possible constants
impl_rational_arithmetics!(DEFAULT_LIMBS, DEFAULT_LIMBS, { 2 * DEFAULT_LIMBS });
impl_rational_arithmetics!(DEFAULT_LIMBS, { 2 * DEFAULT_LIMBS }, { 3 * DEFAULT_LIMBS });
impl_rational_arithmetics!({ 2 * DEFAULT_LIMBS }, { 2 * DEFAULT_LIMBS }, {
    4 * DEFAULT_LIMBS
});

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)   0 <= num   <= n_max,
///  ii)  0 <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
impl<const L: usize> From<HenselCode<L>> for Rational<L> {
    fn from(hc: HenselCode) -> Self {
        let n_max = (&(&hc.modulus() - 1) / 2).sqrt();

        // perform (modified) extended euclidean algorithm on (g, n % g)
        let (mut x0, mut x1) = (hc.modulus(), hc.to_bigint());
        let (mut y0, mut y1) = (BigInt::from(0), BigInt::from(1));
        while x0 > n_max {
            let q = &x0 / &x1;
            let new_x0 = x1.clone();
            let big_x1 = x0 - (&q * &x1);
            let new_y0 = y1.clone();
            let big_y1 = y0 - (&q * &y1);
            (x0, x1) = (new_x0, big_x1);
            (y0, y1) = (new_y0, big_y1);
        }

        Rational::<L> { num: x0, denom: y0 }
    }
}
