use super::{
    fmt,
    ops::{Add, Mul},
    rational::Rational,
};
use crate::bigint::BigInt;
use crate::shared::{Bounded, DEFAULT_LIMBS};

use crypto_bigint::modular::runtime_mod::{DynResidue, DynResidueParams};

// the operation `chinese_remainder` changes the size of the modulus, so we need to track it using a const generics LIMBS
#[derive(Clone)]
pub struct HenselCode<const L: usize = DEFAULT_LIMBS> {
    params: DynResidueParams<L>,
    res: DynResidue<L>, // internal variable that stores the residue
}
impl<const L: usize> HenselCode<L> {
    /// return a BigInt n, with residue `res` mod `params.modulus()`
    pub fn to_bigint(&self) -> BigInt<L> {
        BigInt::new(self.res.retrieve())
    }

    /// return the modulus
    pub fn modulus(&self) -> BigInt<L> {
        BigInt::new(*self.params.modulus())
    }
}

impl<const L: usize> Bounded for HenselCode<L> {
    const L: usize = L;
}

/// Create an HenselCode from two BigInt
pub fn new_hensel_code<const LG: usize, const LN: usize>(
    g: &BigInt<LG>,
    n: &BigInt<LN>,
) -> HenselCode<LG> {
    let params = DynResidueParams::new(&g.0 .0);
    let res = DynResidue::new(&(&n.resize::<LG>() % g).0 .0, params);
    HenselCode { params, res }
}

impl<const L: usize> HenselCode<L> {
    pub fn invert(&self) -> HenselCode<L> {
        HenselCode {
            params: self.params,
            res: self.res.invert().0,
        }
    }
}
/// add two HenselCodes
impl<const L: usize> Add<HenselCode<L>> for HenselCode<L> {
    type Output = HenselCode<L>;
    fn add(self, other: HenselCode<L>) -> HenselCode<L> {
        &self + &other
    }
}
/// multiply two HenselCodes
impl<const L: usize> Mul<HenselCode<L>> for HenselCode<L> {
    type Output = HenselCode<L>;
    fn mul(self, other: HenselCode<L>) -> HenselCode<L> {
        &self * &other
    }
}

/// add two &HenselCodes
impl<'a, 'b, const L: usize> Add<&'b HenselCode<L>> for &'a HenselCode<L> {
    type Output = HenselCode<L>;
    fn add(self, other: &'b HenselCode<L>) -> HenselCode<L> {
        if self.modulus() != other.modulus() {
            panic!("cannot add '{}' and '{}'", self, other);
        }
        HenselCode {
            params: self.params,
            res: self.res + other.res,
        }
    }
}
/// multiply two &HenselCodes
impl<'a, 'b, const L: usize> Mul<&'b HenselCode<L>> for &'a HenselCode<L> {
    type Output = HenselCode<L>;
    fn mul(self, other: &'b HenselCode<L>) -> HenselCode<L> {
        if self.modulus() != other.modulus() {
            panic!("cannot add '{}' and '{}'", self, other);
        }
        HenselCode {
            params: self.params,
            res: self.res * other.res,
        }
    }
}

pub fn chinese_remainder<const L: usize>(hc1: HenselCode<L>, hc2: HenselCode<L>) -> HenselCode<L> {
    let (g1, n1) = (hc1.modulus(), hc1.to_bigint());
    let (g2, n2) = (hc2.modulus(), hc2.to_bigint());
    let g12 = hc1.modulus() * hc2.modulus();
    let (residue_params1, residue_params2, residue_params) = (
        DynResidueParams::new(&g1.0 .0),
        DynResidueParams::new(&g2.0 .0),
        DynResidueParams::new(&g12.0 .0),
    );
    let (mut res_g1, mut res_g2) = (
        DynResidue::new(&g1.0 .0, residue_params2),
        DynResidue::new(&g2.0 .0, residue_params1),
    );
    // i1*g1 = 1 (mod g2), i2*g2 = 1 (mod g1)
    // we need to convert i1 and i2 to a residue mod g1*g2
    let (i1, i2) = (
        DynResidue::new(&res_g1.invert().0.retrieve(), residue_params),
        DynResidue::new(&res_g2.invert().0.retrieve(), residue_params),
    );

    // change modulus g1 -> g1*g2 and g2 -> g1*g2
    (res_g1, res_g2) = (
        DynResidue::new(&g1.0 .0, residue_params),
        DynResidue::new(&g2.0 .0, residue_params),
    );

    let (res_n1, res_n2) = (
        DynResidue::new(&n1.0 .0, residue_params),
        DynResidue::new(&n2.0 .0, residue_params),
    );

    let res = res_g1 * i1 * res_n2 + res_g2 * i2 * res_n1;

    HenselCode {
        params: residue_params,
        res,
    }
}

/// pretty-print HenselCode
impl<const L: usize> fmt::Display for HenselCode<L> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} (mod {})", self.to_bigint(), self.modulus())
    }
}

/// given a prime `p` and a rational `r = num/denom`, where p does not divide denom,
/// returns `r (mod p)`
impl<const L: usize> From<(&BigInt<L>, &Rational<L>)> for HenselCode<L> {
    fn from(params: (&BigInt<L>, &Rational<L>)) -> Self {
        let (g, r) = params;
        let params = DynResidueParams::new(&g.0 .0);
        let denom = DynResidue::<L>::new(&r.denom.0 .0, params);
        let num = DynResidue::<L>::new(&r.num.0 .0, params);
        let (id, _) = denom.invert();
        let res = id * num;
        HenselCode { params, res }
    }
}
