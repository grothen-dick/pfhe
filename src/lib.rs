pub mod bigint;
pub mod crypto_parameters;
pub mod hensel_code;
pub mod macros;
pub mod rational;
pub mod shared;

use std::{clone::Clone, fmt, ops};

#[cfg(test)]
mod tests {
    use crate::crypto_parameters::{EncryptionScheme, PrivateKeySchemeCryptographicParameters};

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
            denom: T::one(),
        };
        simple_tester(&r1, &p);

        // integer inverse
        let r2 = Rational::<T> {
            num: T::one(),
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
            denom: T::one(),
        };
        simple_tester(&r1, &p);

        // integer inverse
        let r2 = Rational::<T> {
            num: T::one(),
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
        let crypto_params: PrivateKeySchemeCryptographicParameters<T> =
            PrivateKeySchemeCryptographicParameters::<T>::new(p1, p2, p3, p4, p5);
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

    // #[test]
    // fn encrypt_add_decrypt() {
    //     let (p1, p2, p3, p4, p5) = (
    //         T::from_u128(7919),
    //         T::from_u128(37),
    //         T::from_u128(41),
    //         T::from_u128(5897),
    //         T::from_u128(7759),
    //     );
    //     let crypto_param = PrivateKeySchemeCryptographicParameters::new(
    //         p1.clone(),
    //         p2.clone(),
    //         p3.clone(),
    //         p4.clone(),
    //         p5.clone(),
    //     );
    //     let r1 = Rational::<T> {
    //         num: T::from_u128(37),
    //         denom: T::from_u128(7),
    //     };
    //     let r2 = Rational::<T> {
    //         num: T::from_u128(43),
    //         denom: T::from_u128(7),
    //     };
    //     let r_sum = (&r1 + &r2).reduce();
    //     println!("rational 1: {r1}");
    //     println!("rational 2: {r2}");
    //     let encrypted_1 = crypto_param.encrypt(r1);
    //     let encrypted_2 = crypto_param.encrypt(r2);
    //     println!("encrypted message 1: {encrypted_1}");
    //     println!("encrypted message 2: {encrypted_2}");
    //     let encrypted_sum = encrypted_1 + encrypted_2;
    //     let sum = crypto_param.decrypt(encrypted_sum);
    //     assert_eq!(sum, r_sum);
    // }
}
