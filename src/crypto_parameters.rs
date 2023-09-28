extern crate crypto_bigint;

use crate::{
    bigint::BigInt,
    hensel_code::{chinese_remainder, new_hensel_code, HenselCode},
    rational::Rational,
    shared::Bounded,
};

use std::convert::From;

/// This is a private key, with five private parameters.
/// Rust doesn't like "const generics expressions" so it is needed to assume that
/// the product p1*...*p5 is representable by a BigInt of size L.
pub struct CryptographicParameters<const L: usize> {
    _p1: BigInt<L>,
    _p2: BigInt<L>,
    _p3: BigInt<L>,
    _p4: BigInt<L>,
    _p5: BigInt<L>,
}

impl<const L: usize> Bounded for CryptographicParameters<L> {
    const L: usize = L;
}

impl<const L: usize> CryptographicParameters<L> {
    pub fn new(
        _p1: BigInt<L>,
        _p2: BigInt<L>,
        _p3: BigInt<L>,
        _p4: BigInt<L>,
        _p5: BigInt<L>,
    ) -> CryptographicParameters<L> {
        CryptographicParameters::<L> {
            _p1,
            _p2,
            _p3,
            _p4,
            _p5,
        }
    }

    /// Returns the product of the 5 primes used as crypto parameters
    pub fn public_key(&self) -> BigInt<L> {
        self._p2 * self._p3 * self._p4 * self._p5
    }

    pub fn chinese_remainder(&self, n1: BigInt<L>, n2: BigInt<L>, n3: BigInt<L>) -> HenselCode<L> {
        let hc1 = new_hensel_code(&self._p1, &n1);
        let hc2 = new_hensel_code(&self._p2, &n2);
        let hc3 = new_hensel_code(&self._p3, &n3);
        chinese_remainder(chinese_remainder(hc1, hc2), hc3)
    }

    pub fn encrypt(&self, m: Rational<L>) -> HenselCode<L> {
        let delta_max: BigInt<L> = self._p1 * self._p2 * self._p3 * self._p5;
        let g: BigInt<L> = delta_max * self._p4;
        let s1 = BigInt::<L>::random_mod(&self._p1);
        let s2 = BigInt::<L>::random_mod(&self._p2);
        let s3 = BigInt::<L>::random_mod(&self._p3);
        let delta = BigInt::<L>::random_mod(&delta_max);

        let dp4: HenselCode<L> = new_hensel_code(&g, &(delta * self._p4));
        let zero = BigInt::<L>::from(0);
        let one = BigInt::<L>::from(1);
        println!("s1: {}", s1);
        println!("s2: {}", s2);
        println!("s3: {}", s3);

        let rs1 = Rational {
            num: s1,
            denom: one,
        };
        // intermediary step with a Farey fraction
        let hc_noise = self.chinese_remainder(zero, s2, s3);
        println!("hc_noise: {}", hc_noise);
        println!(
            "gcd(noise, p1): {}",
            BigInt::gcd(&hc_noise.to_bigint(), &self._p1)
        );
        println!("delta*p4: {}", dp4);
        println!(
            "gcd(delta*p4, p4): {}",
            BigInt::gcd(&dp4.to_bigint(), &self._p4)
        );

        let r_noise: Rational<L> = Rational::<L>::from(&hc_noise);
        println!("r_noise: {}", r_noise);
        let mut rational_term: Rational<L> = rs1 * r_noise;
        println!("rational_term = s1*r_noise: {}", rational_term);
        rational_term = rational_term + m;
        println!("rational_term + m: {}", rational_term);
        // returns a HenselCode
        HenselCode::from((&g, &rational_term)) + dp4
    }

    pub fn decrypt(&self, hc: HenselCode<L>) -> Rational<L> {
        let hc_p4 = new_hensel_code(&self._p4, &hc.to_bigint());
        let r_p4: Rational<L> = Rational::<L>::from(&hc_p4);
        Rational::<L>::from(&HenselCode::<L>::from((&self._p1, &r_p4)))
    }
}

#[cfg(test)]
mod tests {
    use super::CryptographicParameters;
    use crate::hensel_code;

    type BigInt = crate::bigint::BigInt;

    #[test]
    fn chinese_remainder() {
        let (p1, p2, p3, p4, p5) = (
            BigInt::from(4919),
            BigInt::from(7),
            BigInt::from(11),
            BigInt::from(13),
            BigInt::from(17),
        );
        let crypto_param = CryptographicParameters::new(
            p1.clone(),
            p2.clone(),
            p3.clone(),
            p4.clone(),
            p5.clone(),
        );
        let (n1, n2, n3) = (BigInt::from(38), BigInt::from(2), BigInt::from(1));
        let result = crypto_param.chinese_remainder(n1.clone(), n2.clone(), n3.clone());

        assert_eq!((result.to_bigint() % BigInt::from(4919)), n1.clone());
        assert_eq!((result.to_bigint() % BigInt::from(7)), n2.clone());
        assert_eq!((result.to_bigint() % BigInt::from(11)), n3.clone());

        let hc1 = hensel_code::new_hensel_code(&p1, &n1);
        let hc2 = hensel_code::new_hensel_code(&p2, &n2);
        let hc3 = hensel_code::new_hensel_code(&p3, &n3);
        let hc12 = hensel_code::chinese_remainder(hc1, hc2);
        let hc = hensel_code::chinese_remainder(hc12, hc3);
        assert_eq!(result.to_bigint(), hc.to_bigint());
        println!("{} : {}", hc, result);
    }
}
