extern crate crypto_bigint;

use crate::{
    bigint::BigInt,
    hensel_code::{new_hensel_code, HenselCode},
    rational::Rational,
    shared::Bounded,
};

use std::{
    convert::From,
    ops::{Add, Mul},
};

/// This is a private key, with five private parameters
struct CryptographicParameters<const L: usize> {
    p1: BigInt<L>,
    p2: BigInt<L>,
    p3: BigInt<L>,
    p4: BigInt<L>,
    p5: BigInt<L>,
}

impl<const L: usize> Bounded for CryptographicParameters<L> {
    const L: usize = L;
}

// A LOT of trait bounds coming up...
// this is just to tell the compiler that we can indeed do
// arithmetics and hensel_code <-> rational conversions, and tell it to keep track of the
// size L (it's a bit dumb so we need to be quite verbose)
impl<'a, const L: usize> CryptographicParameters<L>
where
    [(); L + L]:,
    [(); L + L + L]:,
    [(); L + L + L + L]:,
    [(); L + L + L + L + L]:,
    BigInt<{ L + L + L + L }>: Sized,
    BigInt<{ L + L + L + L + L }>: Sized,
    Rational<{ L }>: Mul<Rational<{ L + L + L }>, Output = Rational<{ L + L + L + L }>>,
    Rational<{ L + L + L }>: Sized + From<HenselCode<{ L + L + L }>>,
    Rational<{ L + L + L + L }>:
        Sized + Add<Rational<{ L + L + L + L }>, Output = Rational<{ L + L + L + L }>>,
    Rational<{ L + L + L + L + L }>: Sized,
    HenselCode<{ L + L + L }>: Sized,
    HenselCode<{ L + L + L + L + L }>: Sized
        + From<(
            BigInt<{ L + L + L + L + L }>,
            Rational<{ L + L + L + L + L }>,
        )>,
{
    pub fn new(
        p1: BigInt<L>,
        p2: BigInt<L>,
        p3: BigInt<L>,
        p4: BigInt<L>,
        p5: BigInt<L>,
    ) -> CryptographicParameters<L> {
        CryptographicParameters::<L> { p1, p2, p3, p4, p5 }
    }

    /// return the product of the 5 primes used as crypto parameters
    pub fn public_key(&self) -> BigInt<{ L + L + L + L + L }> {
        self.p1.clone() * self.p2.clone() * self.p3.clone() * self.p4.clone() * self.p5.clone()
    }

    pub fn chinese_remainder(
        &self,
        n1: BigInt<{ L }>,
        n2: BigInt<{ L }>,
        n3: BigInt<{ L }>,
    ) -> HenselCode<{ L + L + L }> {
        let hc1 = new_hensel_code(self.p1.clone(), n1);
        let hc2 = new_hensel_code(self.p2.clone(), n2);
        let hc3 = new_hensel_code(self.p3.clone(), n3);
        HenselCode::<{ L + L + L }>::chinese_remainder(
            HenselCode::<{ L + L }>::chinese_remainder(hc1, hc2),
            hc3,
        )
    }

    pub fn encode(&self, m: BigInt<L>) -> HenselCode<{ L + L + L + L + L }> {
        let delta_max: BigInt<{ L + L + L + L }> =
            self.p1.clone() * self.p2.clone() * self.p3.clone() * self.p5.clone();
        let g: BigInt<{ L + L + L + L + L }> = &delta_max * &self.p4;
        let s1 = BigInt::<L>::random_mod(&self.p1);
        let s2 = BigInt::<L>::random_mod(&self.p2);
        let s3 = BigInt::<L>::random_mod(&self.p3);
        let delta = BigInt::<{ L + L + L + L }>::random_mod(&delta_max);

        let dp4: HenselCode<{ L + L + L + L + L }> = new_hensel_code(g.clone(), &delta * &self.p4);
        let zero = BigInt::<L>::from(0);
        let one = BigInt::<L>::from(1);

        let rm = Rational {
            num: m.resize::<{ L + L + L + L }>(),
            denom: one.resize::<{ L + L + L + L }>(),
        };
        let rs1 = Rational {
            num: s1,
            denom: one,
        };
        // intermediary step with a Farey fraction
        let r_noise: Rational<{ L + L + L }> =
            Rational::<{ L + L + L }>::from(self.chinese_remainder(zero, s2, s3));
        let mut rational_term: Rational<{ L + L + L + L }> = rs1 * r_noise;
        rational_term = rational_term + rm;
        let bigger_r_term = rational_term.resize::<{ L + L + L + L + L }>();
        // returns a HenselCode
        return HenselCode::from((g, bigger_r_term)) + dp4;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U256;

    type BigInt = crate::bigint::BigInt;
    const L: usize = U256::LIMBS;

    #[test]
    fn chinese_remainder() {
        let (p1, p2, p3, p4, p5) = (
            BigInt::from(4919),
            BigInt::from(7),
            BigInt::from(11),
            BigInt::from(13),
            BigInt::from(17),
        );
        let crypto_param = CryptographicParameters::new(p1, p2, p3, p4, p5);
        let (n1, n2, n3) = (BigInt::from(38), BigInt::from(2), BigInt::from(1));
        let result = crypto_param.chinese_remainder(n1, n2, n3);

        assert_eq!((result.to_bigint() % BigInt::from(4919)).0 .0, n1.0 .0);
        assert_eq!((result.to_bigint() % BigInt::from(7)).0 .0, n2.0 .0);
        assert_eq!((result.to_bigint() % BigInt::from(11)).0 .0, n3.0 .0);

        let hc1 = new_hensel_code(p1, n1);
        let hc2 = new_hensel_code(p2, n2);
        let hc3 = new_hensel_code(p3, n3);
        let hc12 = HenselCode::<L>::chinese_remainder(hc1, hc2);
        let hc = HenselCode::<L>::chinese_remainder(hc12, hc3);
        assert_eq!(result.to_bigint().0 .0, hc.to_bigint().0 .0);
        println!("{} : {}", hc, result);
    }
}
