extern crate crypto_bigint;

use crate::{
    bigint::BigIntTrait,
    hensel_code::{chinese_remainder, new_hensel_code, HenselCode},
    rational::Rational,
};

use std::convert::From;

/// This is a private key, with five private parameters.
/// Rust doesn't like "const generics expressions" so it is needed to assume that
/// the product p1*...*p5 is representable by a BigInt of size L.
pub struct PrivateKeySchemeCryptographicParameters<T: BigIntTrait> {
    _p1: T,
    _p2: T,
    _p3: T,
    _p4: T,
    _p5: T,
}

pub trait EncryptionScheme<T: BigIntTrait> {
    fn encrypt(&self, m: Rational<T>) -> HenselCode<T>;
    fn decrypt(&self, hc: HenselCode<T>) -> Rational<T>;
}

impl<T: BigIntTrait> PrivateKeySchemeCryptographicParameters<T> {
    pub fn new(_p1: T, _p2: T, _p3: T, _p4: T, _p5: T) -> Self {
        Self {
            _p1,
            _p2,
            _p3,
            _p4,
            _p5,
        }
    }

    /// generates 5 distincts primes from security parameters `lambda, d`
    pub fn new_from_params(lambda: u32, d: u32) -> Self {
        let rho = lambda;
        let eta = 2 * (d + 2) * lambda;
        let gamma: u32 = (lambda / lambda.ilog2()) * (eta - rho).pow(2);
        let mu = gamma - eta - 2 * lambda;
        let mut primes: Vec<T> = Vec::new();
        for size in [(rho + 1), (rho / 2), (rho / 2), eta, mu] {
            loop {
                let current_p = T::generate_prime(Some(size as usize));
                if !primes.contains(&current_p) {
                    primes.push(current_p);
                    break;
                }
            }
        }
        Self::new(
            primes[0].clone(),
            primes[1].clone(),
            primes[2].clone(),
            primes[3].clone(),
            primes[4].clone(),
        )
    }

    /// returns a number `n` such that `n = n1 (mod p1)`, `n = n2 (mod p2)`, `n = n3 (mod p3)`
    fn chinese_remainder(&self, n1: T, n2: T, n3: T) -> HenselCode<T> {
        let hc1 = new_hensel_code(&self._p1, &n1);
        let hc2 = new_hensel_code(&self._p2, &n2);
        let hc3 = new_hensel_code(&self._p3, &n3);
        chinese_remainder(chinese_remainder(hc1, hc2), hc3)
    }
}

impl<T: BigIntTrait> EncryptionScheme<T> for PrivateKeySchemeCryptographicParameters<T> {
    fn encrypt(&self, m: Rational<T>) -> HenselCode<T> {
        let delta_max: T = self._p1.mul(&self._p2).mul(&self._p3).mul(&self._p5);
        let g: T = delta_max.mul(&self._p4);
        let s1 = T::random_mod(&self._p1);
        let s2 = T::random_mod(&self._p2);
        let s3 = T::random_mod(&self._p3);
        let delta = T::random_mod(&delta_max);

        let dp4: HenselCode<T> = new_hensel_code(&g, &(delta.mul(&self._p4)));
        let zero = T::from_u128(0);
        let one = T::from_u128(1);

        // generate an encoding of zero
        let hc_noise = self.chinese_remainder(zero, s2, s3);
        // divide the result by p1 in order to get a correct HenselCode -> Rational conversion
        let hc_noise_1 = HenselCode::<T>::from((
            &(self._p1.mul(&self._p2).mul(&self._p3)),
            &Rational::<T> {
                num: hc_noise.res,
                denom: self._p1.clone(),
            },
        ));

        // convert to a Rational
        let r_noise: Rational<T> = Rational::<T> {
            num: self._p1.clone(),
            denom: T::from_u128(1),
        } * Rational::<T>::from(&hc_noise_1);

        // create a Rational from s1
        let rs1 = Rational {
            num: s1,
            denom: one,
        };
        // multiply rational encoding of zero by s1
        let mut rational_term: Rational<T> = rs1 * r_noise;

        // add the message `m` (a Rational by assumption)
        rational_term = rational_term + m;

        // convert to HenselCode, add another noise `delta*p4`
        // return the result
        HenselCode::from((&g, &rational_term)) + dp4
    }

    fn decrypt(&self, hc: HenselCode<T>) -> Rational<T> {
        let hc_p4 = new_hensel_code(&self._p4, &hc.res);
        let r_p4: Rational<T> = Rational::<T>::from(&hc_p4);
        Rational::<T>::from(&HenselCode::<T>::from((&self._p1, &r_p4)))
    }
}

/// This is a private key, with five private parameters.
/// Rust doesn't like "const generics expressions" so it is needed to assume that
/// the product p1*...*p5 is representable by a BigInt of size L.
pub struct PublicKeySchemeCryptographicParameters<T: BigIntTrait> {
    _p1: T,
    _p2: T,
    _p3: T,
    _p4: T,
    lambda: u32,
    g: T,
    g_prime: T,
    e: HenselCode<T>,
}

impl<T: BigIntTrait> PublicKeySchemeCryptographicParameters<T> {
    pub fn new(_p1: T, _p2: T, _p3: T, _p4: T, lambda: u32, e: HenselCode<T>) -> Self {
        let g = _p1.mul(&_p2.mul(&_p3.mul(&_p4)));
        let g_prime = _p3.mul(&_p4);
        Self {
            _p1,
            _p2,
            _p3,
            _p4,
            lambda,
            g,
            g_prime,
            e,
        }
    }

    /// generates 5 distincts primes from security parameters `lambda, d`
    pub fn new_from_params(lambda: u32, d: u32) -> Self {
        let rho = lambda;
        let eta = d * lambda;
        let mu = d.pow(2) * lambda * lambda.ilog(2) - eta - 2 * lambda - 3;
        let gamma = 2 * eta + 3 * lambda / 2 + mu + 3;
        let mut primes: Vec<T> = Vec::new();
        for size in [rho, rho + 3, eta, eta] {
            loop {
                let current_p = T::generate_prime(Some(size as usize));
                if !primes.contains(&current_p) {
                    primes.push(current_p);
                    break;
                }
            }
        }
        let [p1, p2, p3, p4prime] = &primes[..] else {
            todo!()
        };
        let p4 = p4prime.pow((f64::from(mu) / f64::from(eta + 1)).ceil() as u128);
        let g = p1.mul(&p2.mul(&p3.mul(&p4)));
        let t = T::random_mod(&T::from_u128(2u128.pow(lambda - 1)));
        let delta_e = T::random_mod(&T::from_u128(2u128.pow(gamma - eta)));
        let hc_noise = chinese_remainder(
            new_hensel_code(p1, &T::from_u128(0)),
            new_hensel_code(p2, &t),
        );
        let hc_p3_res = HenselCode::<T>::from((p3, &Rational::from(&hc_noise))).res;
        let e = new_hensel_code(&g, &hc_p3_res.add(&delta_e.mul(p3)));
        Self::new(
            primes[0].clone(),
            primes[1].clone(),
            primes[2].clone(),
            primes[3].clone(),
            lambda,
            e,
        )
    }
}

impl<T: BigIntTrait> EncryptionScheme<T> for PublicKeySchemeCryptographicParameters<T> {
    fn encrypt(&self, m: Rational<T>) -> HenselCode<T> {
        let s1: T = T::random_mod(&T::from_u128(2).pow((self.lambda - 1) as u128));
        let s2: T = T::random_mod(&T::from_u128(2).pow((self.lambda - 1) as u128));
        let p12 = self._p1.mul(&self._p2);
        let delta: T = T::random_mod(&p12.mul(&self._p4).sub(&p12)).add(&p12);
        let hc_prime_res = (new_hensel_code(&self.g_prime, &s1.mul(&self.e.res))
            + HenselCode::from((&self.g_prime, &m)))
        .res;
        let encrypted_res = hc_prime_res
            .add(&s2.mul(&self.g_prime))
            .add(&delta.mul(&self.g.pow(2)));
        new_hensel_code(&self.g, &encrypted_res)
    }

    fn decrypt(&self, hc: HenselCode<T>) -> Rational<T> {
        let hc_p3 = new_hensel_code(&self._p3, &hc.res);
        let r_p3: Rational<T> = Rational::<T>::from(&hc_p3);
        Rational::<T>::from(&HenselCode::<T>::from((&self._p1, &r_p3)))
    }
}

#[cfg(test)]
mod tests {
    use super::PrivateKeySchemeCryptographicParameters;
    use crate::bigint::BigIntTrait;
    use crate::hensel_code;

    type T = crate::bigint::WrappingCryptoBigInt;

    #[test]
    fn chinese_remainder() {
        let _from_u128 = <T as BigIntTrait>::from_u128;
        let (p1, p2, p3, p4, p5) = (
            T::from_u128(4919),
            T::from_u128(7),
            T::from_u128(11),
            T::from_u128(13),
            T::from_u128(17),
        );
        let crypto_param = PrivateKeySchemeCryptographicParameters {
            _p1: p1.clone(),
            _p2: p2.clone(),
            _p3: p3.clone(),
            _p4: p4.clone(),
            _p5: p5.clone(),
        };
        let (n1, n2, n3) = (T::from_u128(38), T::from_u128(2), T::from_u128(1));
        let result = crypto_param.chinese_remainder(n1.clone(), n2.clone(), n3.clone());

        assert_eq!((result.res.rem(&T::from_u128(4919))), n1.clone());
        assert_eq!((result.res.rem(&T::from_u128(7))), n2.clone());
        assert_eq!((result.res.rem(&T::from_u128(11))), n3.clone());

        let hc1 = hensel_code::new_hensel_code(&p1, &n1);
        let hc2 = hensel_code::new_hensel_code(&p2, &n2);
        let hc3 = hensel_code::new_hensel_code(&p3, &n3);
        let hc12 = hensel_code::chinese_remainder(hc1, hc2);
        let hc = hensel_code::chinese_remainder(hc12, hc3);
        assert_eq!(result.res, hc.res);
        println!("{} : {}", hc, result);
    }
}
