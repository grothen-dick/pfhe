extern crate crypto_bigint;

use crate::{
    big_int::BigInt,
    hensel_code::{new_hensel_code, HenselCode},
    rational::Rational,
};

/// This is a private key, with five private parameters
struct CryptographicParameters {
    p1: BigInt,
    p2: BigInt,
    p3: BigInt,
    p4: BigInt,
    p5: BigInt,
}

impl CryptographicParameters {
    pub fn new(
        p1: BigInt,
        p2: BigInt,
        p3: BigInt,
        p4: BigInt,
        p5: BigInt,
    ) -> CryptographicParameters {
        CryptographicParameters { p1, p2, p3, p4, p5 }
    }
    /// return the product of the 5 primes used as crypto parameters
    pub fn public_key(&self) -> BigInt {
        self.p1 * self.p2 * self.p3 * self.p4 * self.p5
    }

    pub fn chinese_remainder(&self, n1: BigInt, n2: BigInt, n3: BigInt) -> HenselCode {
        let hc1 = new_hensel_code(self.p1, n1);
        let hc2 = new_hensel_code(self.p2, n2);
        let hc3 = new_hensel_code(self.p3, n3);
        HenselCode::chinese_remainder(HenselCode::chinese_remainder(hc1, hc2), hc3)
    }

    pub fn encode(&self, m: BigInt) -> HenselCode {
        let delta_max = self.p1 * self.p2 * self.p3 * self.p5;
        let g = delta_max * self.p4;
        let s1 = BigInt::random_mod(self.p1);
        let s2 = BigInt::random_mod(self.p2);
        let s3 = BigInt::random_mod(self.p3);
        let delta = BigInt::random_mod(delta_max);

        let dp4 = new_hensel_code(g, delta * self.p4);

        let rm = Rational { num: m, denom: 1 };
        let rs1 = Rational { num: s1, denom: 1 };
        // intermediary step with a Farey fraction
        let rational_term = s1 * Rational::from(self.chinese_remainder(0, s2, s3)) + rm;
        // returns a HenselCode
        return HenselCode::from((g, rational_term)) + dp4;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::U256;
    type BigInt = crate::big_int::BigInt;
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
        // for 0 reason the usual 'modulo' operator `%` has a range of [-p, p]
        // for `(...) % p` instead of the normal range [0, p]. The sensible version
        // is instead written `(...).rem_euclid(p)`.
        assert_eq!(result.n.rem(U256::from(4919)), n1);
        assert_eq!(result.n.rem(U256::from(7)), n2);
        assert_eq!(result.n.rem(U256::from(11)), n3);
        let hc1 = new_hensel_code(p1, n1);
        let hc2 = new_hensel_code(p2, n2);
        let hc3 = new_hensel_code(p3, n3);
        let hc12 = HenselCode::<L>::chinese_remainder(hc1, hc2);
        let hc = HenselCode::<L>::chinese_remainder(hc12, hc3);
        assert_eq!(result.n.rem_euclid(p1 * p2 * p3), hc.n);
        println!("{} : {}", hc, result);
    }
}
