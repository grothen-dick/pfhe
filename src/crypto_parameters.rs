extern crate crypto_bigint;

use crate::{
    bigint::BigInt,
    hensel_code::{chinese_remainder, new_hensel_code, HenselCode},
    rational::Rational,
    shared::Bounded,
};

use crypto_primes::generate_prime;
use std::convert::From;

pub trait EncryptionScheme<const L: usize> {
    fn chinese_remainder(&self, n1: BigInt<L>, n2: BigInt<L>, n3: BigInt<L>) -> HenselCode<L>;
    fn encrypt(&self, m: Rational<L>) -> HenselCode<L>;
    fn decrypt(&self, hc: HenselCode<L>) -> Rational<L>;
}

/// This is a private key, with five private parameters.
/// Rust doesn't like "const generics expressions" so it is needed to assume that
/// the product p1*...*p5 is representable by a BigInt of size L.
pub struct PrivateKeySchemeCryptographicParameters<const L: usize> {
    _p1: BigInt<L>,
    _p2: BigInt<L>,
    _p3: BigInt<L>,
    _p4: BigInt<L>,
    _p5: BigInt<L>,
}

impl<const L: usize> PrivateKeySchemeCryptographicParameters<L> {
    pub fn new(lambda: u32, d: u32) -> PrivateKeySchemeCryptographicParameters<L> {
        let rho = lambda;
        let eta = 2 * (d + 2) * lambda;
        let gamma: u32 = (lambda / lambda.ilog2()) * (eta - rho).pow(2);
        let mu = gamma - eta - 2 * lambda;
        PrivateKeySchemeCryptographicParameters::<L> {
            _p1: BigInt::new(generate_prime(Some(((rho + 1) as u16).into()))),
            _p2: BigInt::new(generate_prime(Some(((rho / 2) as u16).into()))),
            _p3: BigInt::new(generate_prime(Some(((rho / 2) as u16).into()))),
            _p4: BigInt::new(generate_prime(Some((eta as u16).into()))),
            _p5: BigInt::new(generate_prime(Some((mu as u16).into()))),
        }
    }
}

impl<const L: usize> Bounded for PrivateKeySchemeCryptographicParameters<L> {
    const L: usize = L;
}

impl<const L: usize> EncryptionScheme<L> for PrivateKeySchemeCryptographicParameters<L> {
    /// returns a number `n` such that `n = n1 (mod p1)`, `n = n2 (mod p2)`, `n = n3 (mod p3)`
    fn chinese_remainder(&self, n1: BigInt<L>, n2: BigInt<L>, n3: BigInt<L>) -> HenselCode<L> {
        let hc1 = new_hensel_code(&self._p1, &n1);
        let hc2 = new_hensel_code(&self._p2, &n2);
        let hc3 = new_hensel_code(&self._p3, &n3);
        chinese_remainder(chinese_remainder(hc1, hc2), hc3)
    }

    fn encrypt(&self, m: Rational<L>) -> HenselCode<L> {
        let delta_max: BigInt<L> = self._p1 * self._p2 * self._p3 * self._p5;
        let g: BigInt<L> = delta_max * self._p4;
        let s1 = BigInt::<L>::random_mod(&self._p1);
        let s2 = BigInt::<L>::random_mod(&self._p2);
        let s3 = BigInt::<L>::random_mod(&self._p3);
        let delta = BigInt::<L>::random_mod(&delta_max);

        let dp4: HenselCode<L> = new_hensel_code(&g, &(delta * self._p4));
        let zero = BigInt::<L>::from(0);
        let one = BigInt::<L>::from(1);

        // generate an encoding of zero
        let hc_noise = self.chinese_remainder(zero, s2, s3);
        // divide the result by p1 in order to get a correct HenselCode -> Rational conversion
        let hc_noise_1 = HenselCode::<L>::from((
            &(self._p1 * self._p2 * self._p3),
            &Rational::<L> {
                num: hc_noise.to_bigint(),
                denom: self._p1,
            },
        ));

        // convert to a Rational
        let r_noise: Rational<L> = Rational::<L> {
            num: self._p1,
            denom: BigInt::<L>::from(1),
        } * Rational::<L>::from(&hc_noise_1);

        // create a Rational from s1
        let rs1 = Rational {
            num: s1,
            denom: one,
        };
        // multiply rational encoding of zero by s1
        let mut rational_term: Rational<L> = rs1 * r_noise;

        // add the message `m` (a Rational by assumption)
        rational_term = rational_term + m;

        // convert to HenselCode, add another noise `delta*p4`
        // return the result
        HenselCode::from((&g, &rational_term)) + dp4
    }

    fn decrypt(&self, hc: HenselCode<L>) -> Rational<L> {
        let hc_p4 = new_hensel_code(&self._p4, &hc.to_bigint());
        let r_p4: Rational<L> = Rational::<L>::from(&hc_p4);
        Rational::<L>::from(&HenselCode::<L>::from((&self._p1, &r_p4)))
    }
}

#[cfg(test)]
mod tests {
    use super::{EncryptionScheme, PrivateKeySchemeCryptographicParameters};
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
        let crypto_param = PrivateKeySchemeCryptographicParameters {
            _p1: p1.clone(),
            _p2: p2.clone(),
            _p3: p3.clone(),
            _p4: p4.clone(),
            _p5: p5.clone(),
        };
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
