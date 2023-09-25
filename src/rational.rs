use super::{
    fmt,
    hensel_code::HenselCode,
    ops::{Add, Mul},
    Clone,
};

use crate::bigint::BigInt;
use crate::shared::{Bounded, DEFAULT_LIMBS};

pub struct Rational<const L: usize = DEFAULT_LIMBS> {
    pub num: BigInt<L>,
    pub denom: BigInt<L>,
}

impl<const L: usize> Bounded for Rational<L> {
    const L: usize = L;
}

impl<const L: usize> Clone for Rational<L> {
    fn clone(&self) -> Self {
        Rational {
            num: self.num.clone(),
            denom: self.denom.clone(),
        }
    }
}

/// we want to simplify a rational
impl<const L: usize> Rational<L> {
    pub fn reduce(&self) -> Rational<L> {
        let gcd = BigInt::gcd(&self.num, &self.denom);
        let num = &self.num / &gcd;
        let denom = &self.denom / &gcd;
        Rational { num, denom }
    }
    pub fn resize<const Lnew: usize>(&self) -> Rational<Lnew> {
        Rational {
            num: self.num.resize::<Lnew>(),
            denom: self.denom.resize::<Lnew>(),
        }
    }
}

/// Add two &Rational
impl<'a, 'b, const L: usize> Add<&'b Rational<L>> for &'a Rational<L>
where
    [(); L + L]:,
    BigInt<{ L + L }>: Sized,
    Rational<{ L + L }>: Sized,
{
    type Output = Rational<{ L + L }>;
    fn add(self, other: &'b Rational<L>) -> Rational<{ L + L }> {
        let (r1, r2) = (self, other);
        let num1 = &r1.num * &r2.denom;
        let num2 = &r1.denom * &r2.num; // mystical shit: for Rust, L + L != L + L
        let num = num1 + num2;
        let denom: BigInt<{ L + L }> = &r1.denom * &r2.denom;
        Rational { num, denom }.reduce()
    }
}

/// Multiply two &Rational
impl<'a, 'b, const L1: usize, const L2: usize> Mul<&'b Rational<L2>> for &'a Rational<L1>
where
    Rational<{ L1 + L2 }>: Sized,
{
    type Output = Rational<{ L1 + L2 }>;
    fn mul(self, other: &'b Rational<L2>) -> Rational<{ L1 + L2 }> {
        let (r1, r2) = (self, other);
        Rational {
            num: &r1.num * &r2.num,
            denom: &r1.denom * &r2.denom,
        }
        .reduce()
    }
}

/// Add two Rational
impl<const L: usize> Add<Rational<L>> for Rational<L>
where
    [(); L + L]:,
    Rational<{ L + L }>: Sized,
{
    type Output = Rational<{ L + L }>;
    fn add(self, other: Rational<L>) -> Rational<{ L + L }> {
        &self + &other
    }
}

/// Multiply two Rational
impl<const L1: usize, const L2: usize> Mul<Rational<L2>> for Rational<L1>
where
    [(); L1 + L2]:,
    Rational<{ L1 + L2 }>: Sized,
{
    type Output = Rational<{ L1 + L2 }>;
    fn mul(self, other: Rational<L2>) -> Rational<{ L1 + L2 }> {
        let (r1, r2) = (self, other);
        Rational {
            num: &r1.num * &r2.num,
            denom: &r1.denom * &r2.denom,
        }
        .reduce()
    }
}

/// pretty-print Rational
impl<const L: usize> fmt::Display for Rational<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)   0 <= num   <= n_max,
///  ii)  0 <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
impl<const L: usize> From<HenselCode<L>> for Rational<L>
where
    [(); L]:,
    [(); L + L]:,
    BigInt<{ L + L }>: Sized,
{
    fn from(hc: HenselCode<L>) -> Self {
        let n_max = (&(&hc.modulus() - 1) / 2).sqrt();

        // perform (modified) extended euclidean algorithm on (g, n % g)
        let (mut x0, mut x1) = (hc.modulus(), hc.to_bigint());
        let (mut y0, mut y1) = (BigInt::<L>::from(0), BigInt::<L>::from(1));
        while x0 > n_max {
            let q = (&x0 / &x1).resize::<L>();
            let new_x0 = x1.clone();
            let big_x1 = x0.resize::<{ L + L }>() - (&q * &x1);
            let new_y0 = y1.clone();
            let big_y1 = y0.resize::<{ L + L }>() - (&q * &y1);
            (x0, x1) = (new_x0, big_x1.resize::<L>());
            (y0, y1) = (new_y0, big_y1.resize::<L>());
        }

        Rational { num: x0, denom: y0 }
    }
}
