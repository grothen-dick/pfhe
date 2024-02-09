use super::{
    fmt,
    ops::{Add, Mul},
    rational::Rational,
};
use crate::bigint::BigIntTrait;

#[derive(Clone, Debug)]
pub struct HenselCode<T: BigIntTrait> {
    pub modulus: T,
    pub res: T, // internal variable that stores the residue
}
impl<T: BigIntTrait> HenselCode<T> {
    pub fn generate_zero(modulus: T) -> HenselCode<T> {
        HenselCode {
            modulus,
            res: T::from_u128(0),
        }
    }
}

/// Creates an HenselCode from two BigInt
pub fn new_hensel_code<T: BigIntTrait>(modulus: &T, n: &T) -> HenselCode<T> {
    HenselCode {
        modulus: modulus.clone(),
        res: n.rem(modulus),
    }
}

impl<T: BigIntTrait> HenselCode<T> {
    pub fn invert(&self) -> HenselCode<T> {
        let g = self.modulus.clone();
        let (mut x0, mut x1) = (self.modulus.clone(), self.res.clone());
        let (mut z0, mut z1) = (T::from_u128(0), T::from_u128(1));
        while x1.is_zero().unwrap_u8() == 0 {
            let integer_div = x0.div(&x1);
            (x0, x1) = (x1.clone(), x0.sub(&integer_div.mul(&x1)));
            (z0, z1) = (
                z1.clone(),
                (z0.add(&g).sub(&integer_div.mul(&z1).rem(&g))).rem(&g),
            );
        }
        // we have:
        // x0 = gcd(modulus, res) = z0 * res % modulus
        // x1 = 0
        if x0 != T::from_u128(1) {
            HenselCode {
                modulus: self.modulus.clone(),
                res: T::from_u128(0),
            }
        } else {
            HenselCode {
                modulus: self.modulus.clone(),
                res: z0,
            }
        }
    }
}
/// Adds two HenselCodes
impl<T: BigIntTrait> Add<HenselCode<T>> for HenselCode<T> {
    type Output = HenselCode<T>;
    fn add(self, other: HenselCode<T>) -> HenselCode<T> {
        &self + &other
    }
}
/// Multiplies two HenselCodes
impl<T: BigIntTrait> Mul<HenselCode<T>> for HenselCode<T> {
    type Output = HenselCode<T>;
    fn mul(self, other: HenselCode<T>) -> HenselCode<T> {
        &self * &other
    }
}

/// Adds two &HenselCodes
impl<'a, 'b, T: BigIntTrait> Add<&'b HenselCode<T>> for &'a HenselCode<T> {
    type Output = HenselCode<T>;
    fn add(self, other: &'b HenselCode<T>) -> HenselCode<T> {
        if self.modulus != other.modulus {
            panic!("cannot add '{}' and '{}'", self, other);
        }
        HenselCode {
            modulus: self.modulus.clone(),
            res: self.res.add(&other.res).rem(&self.modulus),
        }
    }
}
/// Multiplies two &HenselCodes
impl<'a, 'b, T: BigIntTrait> Mul<&'b HenselCode<T>> for &'a HenselCode<T> {
    type Output = HenselCode<T>;
    fn mul(self, other: &'b HenselCode<T>) -> HenselCode<T> {
        if self.modulus != other.modulus {
            panic!("cannot add '{}' and '{}'", self, other);
        }
        HenselCode {
            modulus: self.modulus.clone(),
            res: self.res.mul(&other.res).rem(&self.modulus),
        }
    }
}

pub fn chinese_remainder<T: BigIntTrait>(hc1: HenselCode<T>, hc2: HenselCode<T>) -> HenselCode<T> {
    let (g1, n1) = (hc1.modulus, hc1.res);
    let (g2, n2) = (hc2.modulus, hc2.res);
    // println!("\ng1: {g1}, g2: {g2}, gcd: {}", g1.gcd(&g2));
    assert!(PartialEq::eq(&g1.gcd(&g2), &T::from_u128(1)));
    let g12 = g1.mul(&g2);
    let (g1_mod_g2, g2_mod_g1) = (
        HenselCode {
            modulus: g2.clone(),
            res: g1.clone(),
        },
        HenselCode {
            modulus: g1.clone(),
            res: g2.clone(),
        },
    );
    // i1*g1 = 1 (mod g2), i2*g2 = 1 (mod g1)
    let (i1, i2) = (g1_mod_g2.invert().res, g2_mod_g1.invert().res);

    new_hensel_code(&g12, &g1.mul(&i1).mul(&n2).add(&g2.mul(&i2).mul(&n1)))
}

/// Pretty-prints HenselCode
impl<T: BigIntTrait> fmt::Display for HenselCode<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} (mod {})", self.res, self.modulus)
    }
}

/// Given a prime `p` and a rational `r = num/denom`, where p does not divide denom,
/// returns `r (mod p)`
impl<T: BigIntTrait> From<(&T, &Rational<T>)> for HenselCode<T> {
    fn from(params: (&T, &Rational<T>)) -> Self {
        let (g, r) = params;

        let denom = new_hensel_code(g, &r.denom);
        let num = new_hensel_code(g, &r.num);

        if PartialEq::ne(&g.gcd(&r.denom), &T::from_u128(1)) {
            Self::generate_zero(g.clone())
        } else {
            num * (denom.invert())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bigint::BigIntTrait;

    type T = crate::bigint::WrappingCryptoBigInt;

    #[test]
    fn invert_hensel_code() {
        let _from_u128 = <T as BigIntTrait>::from_u128;
        let (p1, _p2, _p3) = (T::from_u128(4919), T::from_u128(7), T::from_u128(11));
        let (n1, _n2, _n3) = (T::from_u128(38), T::from_u128(2), T::from_u128(1));

        let hc1 = new_hensel_code(&p1, &n1);

        assert_eq!(hc1.res.mul(&hc1.invert().res).rem(&p1), T::from_u128(1));
    }

    #[test]
    fn chinese_remainder() {
        let (p1, p2) = (T::from_u128(4919), T::from_u128(7));
        let (n1, n2) = (
            new_hensel_code(&p1, &T::from_u128(38)),
            new_hensel_code(&p2, &T::from_u128(2)),
        );
        let result = super::chinese_remainder(n1.clone(), n2.clone());
        assert_eq!((result.res.rem(&T::from_u128(4919))), n1.res.clone());
        assert_eq!((result.res.rem(&T::from_u128(7))), n2.res.clone());
    }
}
