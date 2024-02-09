pub mod bigint;
pub mod crypto_parameters;
pub mod hensel_code;
pub mod macros;
pub mod rational;
pub mod shared;

use std::{clone::Clone, fmt, ops};

#[cfg(test)]
mod tests {
    use crate::crypto_parameters::CryptographicParameters;

    use super::bigint::{BigIntTrait, WrappingCryptoBigInt};
    use super::hensel_code::{new_hensel_code, HenselCode};
    use super::rational::Rational;
    use super::shared::DEFAULT_LIMBS;

    const L: usize = DEFAULT_LIMBS;
    type T = WrappingCryptoBigInt<L>;

    #[test]
    fn translates_rational_to_hensel_code() {
        fn simple_tester(r: &Rational<T>, p: &T) -> () {
            let hc = HenselCode::from((p, r));
            let id_hc = &new_hensel_code(p, &(r.denom)).invert();
            let n_hc = new_hensel_code(p, &(r.num));
            assert_eq!(&hc.modulus, p);
            assert_eq!(hc.res, (id_hc * &n_hc).res);
            println!("rational: {} => hensel code: {}", r, hc);
        }

        let p = T::from_u128(37);

        // positive integer
        let r1 = Rational::<T> {
            num: T::from_u128(6),
            denom: T::from_u128(1),
        };
        simple_tester(&r1, &p);

        // integer inverse
        let r2 = Rational::<T> {
            num: T::from_u128(1),
            denom: T::from_u128(8),
        };
        simple_tester(&r2, &p);

        // general rational
        let r3 = Rational::<T> {
            num: T::from_u128(6),
            denom: T::from_u128(8),
        };
        simple_tester(&r3, &p);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        fn simple_tester(r: &Rational<T>, p: &T) -> () {
            let hc = HenselCode::from((p, r));
            let new_r = Rational::<T>::from(&hc);
            let id_hc = new_hensel_code(p, &r.denom).invert();
            let n_hc = new_hensel_code(p, &r.num);
            assert_eq!(&hc.modulus, p);
            assert_eq!(hc.res, (id_hc * n_hc).res);
            println!(
                "rational: {} => hensel code: {} => rational: {}",
                r, hc, new_r
            );
        }

        let p: T = T::from_u128(7919); // thanks wikipedia for this prime

        // positive integer
        let r1 = Rational::<T> {
            num: T::from_u128(6),
            denom: T::from_u128(1),
        };
        simple_tester(&r1, &p);

        // integer inverse
        let r2 = Rational::<T> {
            num: T::from_u128(1),
            denom: T::from_u128(8),
        };
        simple_tester(&r2, &p);

        // general rational
        let r3 = Rational::<T> {
            num: T::from_u128(6),
            denom: T::from_u128(8),
        };
        simple_tester(&r3, &p);
    }

    #[test]
    fn encrypt_decrypt() {
        let (p1, p2, p3, p4, p5) = (
            T::from_u128(7919),
            T::from_u128(37),
            T::from_u128(41),
            T::from_u128(5897),
            T::from_u128(7759),
        );
        let crypto_params: CryptographicParameters<T> =
            CryptographicParameters::<T>::new(p1, p2, p3, p4, p5);
        let message: Rational<T> = Rational {
            num: T::from_u128(43),
            denom: T::from_u128(44),
        };
        println!("message: {}", message);
        let ciphertext = crypto_params.encrypt(message.clone());
        println!("ciphertext: {}", ciphertext);
        let decrypted = crypto_params.decrypt(ciphertext);
        assert_eq!(message, decrypted);
    }
}
