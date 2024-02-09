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
pub struct CryptographicParameters<T: BigIntTrait> {
    _p1: T,
    _p2: T,
    _p3: T,
    _p4: T,
    _p5: T,
}

impl<T: BigIntTrait> CryptographicParameters<T> {
    pub fn new(_p1: T, _p2: T, _p3: T, _p4: T, _p5: T) -> CryptographicParameters<T> {
        CryptographicParameters::<T> {
            _p1,
            _p2,
            _p3,
            _p4,
            _p5,
        }
    }

    /// generates 5 distincts primes from security parameters `lambda, d`
    pub fn from_params(lambda: u32, d: u32) -> CryptographicParameters<T> {
        let rho = lambda;
        let eta = 2 * (d + 2) * lambda;
        let gamma: u32 = (lambda / lambda.ilog2()) * (eta - rho).pow(2);
        let mu = gamma - eta - 2 * lambda;
        let mut primes: Vec<T> = Vec::new();
        for size in [(rho + 1), (rho / 2), (rho / 2), eta, mu] {
            loop {
                let current_p = T::generate_prime(Some(size as usize));
                // println!("{current_p}");
                // for p in &primes {
                //     print!("{p}, ");
                // }
                // println!("");
                if !primes.contains(&current_p) {
                    primes.push(current_p);
                    break;
                }
            }
        }
        CryptographicParameters::<T> {
            _p1: primes[0].clone(),
            _p2: primes[1].clone(),
            _p3: primes[2].clone(),
            _p4: primes[3].clone(),
            _p5: primes[4].clone(),
        }
    }

    /// Returns the product of the 5 primes used as crypto parameters
    pub fn public_key(&self) -> T {
        self._p1
            .mul(&self._p2)
            .mul(&self._p3)
            .mul(&self._p4)
            .mul(&self._p5)
    }

    /// returns a number `n` such that `n = n1 (mod p1)`, `n = n2 (mod p2)`, `n = n3 (mod p3)`
    pub fn chinese_remainder(&self, n1: T, n2: T, n3: T) -> HenselCode<T> {
        let hc1 = new_hensel_code(&self._p1, &n1);
        let hc2 = new_hensel_code(&self._p2, &n2);
        let hc3 = new_hensel_code(&self._p3, &n3);
        chinese_remainder(chinese_remainder(hc1, hc2), hc3)
    }

    pub fn encrypt(&self, m: Rational<T>) -> HenselCode<T> {
        // println!(
        //     "p1: {}, p2: {}, p3: {}, p4: {}, p5: {}",
        //     self._p1, self._p2, self._p3, self._p4, self._p5
        // );
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

    pub fn decrypt(&self, hc: HenselCode<T>) -> Rational<T> {
        let hc_p4 = new_hensel_code(&self._p4, &hc.res);
        let r_p4: Rational<T> = Rational::<T>::from(&hc_p4);
        Rational::<T>::from(&HenselCode::<T>::from((&self._p1, &r_p4)))
    }
}

#[cfg(test)]
mod tests {
    use super::CryptographicParameters;
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
        let crypto_param = CryptographicParameters::new(
            p1.clone(),
            p2.clone(),
            p3.clone(),
            p4.clone(),
            p5.clone(),
        );
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
