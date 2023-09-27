use super::{
    fmt,
    hensel_code::HenselCode,
    ops::{Add, Mul},
    Clone,
};

use crate::bigint::BigInt;

#[derive(Clone)]
pub struct Rational<const L: usize> {
    pub num: BigInt<L>,
    pub denom: BigInt<L>,
}

/// basic functionalities: simplify common factors, resize
impl<const L: usize> Rational<L> {
    pub fn reduce(&self) -> Self {
        let gcd = BigInt::gcd(&self.num, &self.denom);
        let num = self.num / gcd;
        let denom = self.denom / gcd;
        Rational::<L> { num, denom }
    }
    pub fn resize<const LNEW: usize>(&self) -> Rational<LNEW> {
        Rational::<LNEW> {
            num: self.num.resize::<LNEW>(),
            denom: self.denom.resize::<LNEW>(),
        }
    }
}

/// pretty-print Rational
impl<const L: usize> fmt::Display for Rational<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// Add two &Rational
impl<'a, 'b, const L: usize> Add<&'b Rational<L>> for &'a Rational<L> {
    type Output = Rational<L>;
    fn add(self, other: &'b Rational<L>) -> Rational<L> {
        let (r1, r2) = (self, other);
        let num1 = r1.num * r2.denom;
        let num2 = r1.denom * r2.num;
        let num: BigInt<L> = num1 + num2;
        let denom: BigInt<L> = r1.denom * r2.denom;
        Rational::<L> { num, denom }.reduce()
    }
}

/// Add two Rational
impl<const L: usize> Add<Rational<L>> for Rational<L> {
    type Output = Rational<L>;
    fn add(self, other: Rational<L>) -> Rational<L> {
        &self + &other
    }
}

/// Multiply two &Rational
impl<'a, 'b, const L: usize> Mul<&'b Rational<L>> for &'a Rational<L> {
    type Output = Rational<L>;
    fn mul(self, other: &'b Rational<L>) -> Rational<L> {
        let (r1, r2) = (self, other);
        Rational::<L> {
            num: r1.num * r2.num,
            denom: r1.denom * r2.denom,
        }
        .reduce()
    }
}

/// Multiply two Rational
impl<const L: usize> Mul<Rational<L>> for Rational<L> {
    type Output = Rational<L>;
    fn mul(self, other: Rational<L>) -> Rational<L> {
        &self * &other
    }
}

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)   0 <= num   <= n_max,
///  ii)  0 <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
impl<const L: usize> From<&HenselCode<L>> for Rational<L> {
    fn from(hc: &HenselCode<L>) -> Self {
        let n_max = ((hc.modulus() - BigInt::<L>::from(1)) / BigInt::<L>::from(2)).sqrt();

        // perform (modified) extended euclidean algorithm on (g, n % g)
        let (mut x0, mut x1) = (hc.modulus(), hc.to_bigint());
        let (mut y0, mut y1) = (BigInt::from(0), BigInt::from(1));
        while x0 > n_max {
            let q = x0 / x1;
            let new_x0 = x1;
            let big_x1 = x0 - (q * x1);
            let new_y0 = y1;
            let big_y1 = y0 - (q * y1);
            (x0, x1) = (new_x0, big_x1);
            (y0, y1) = (new_y0, big_y1);
        }

        Rational::<L> { num: x0, denom: y0 }
    }
}
