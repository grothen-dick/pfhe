extern crate lazy_static;
pub mod bigint;
pub mod crypto_parameters;
pub mod hensel_code;
pub mod macros;
pub mod rational;
pub mod shared;

use std::{clone::Clone, fmt, ops};

#[cfg(test)]
mod tests {
    use super::bigint::BigIntTrait;
    use super::hensel_code::{new_hensel_code, HenselCode};
    use super::rational::Rational;
    use crate::crypto_parameters::{
        EncryptionScheme, PrivateKeySchemeCryptographicParameters,
        PublicKeySchemeCryptographicParameters,
    };
    use lazy_static::lazy_static;

    use num_bigint_dig::BigInt;

    //const L: usize = DEFAULT_LIMBS;
    type T = BigInt;

    lazy_static! {
        pub static ref PUBLIC_PARAMS : PublicKeySchemeCryptographicParameters::<T> = PublicKeySchemeCryptographicParameters::<T>::new_from_params(128, 10); // secure parameters would be 512 and 10;
    }

    #[test]
    fn translates_rational_to_hensel_code() {
        fn simple_tester(r: &Rational<T>, p: &T) {
            let hc = HenselCode::from((p, r));
            let id_hc = &new_hensel_code(p, &(r.denom)).invert();
            let n_hc = new_hensel_code(p, &(r.num));
            assert_eq!(&hc.modulus, p);
            assert_eq!(hc.res, (id_hc * &n_hc).res);
            println!("rational: {} => hensel code: {}", r, hc);
        }

        let p = T::from_u128(7919);

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

        let r4 = Rational::<T> {
            num: T::from_u128(32),
            denom: T::from_u128(27),
        };
        simple_tester(&r4, &p);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        fn simple_tester(r: &Rational<T>, p: &T) {
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
            num: T::from_u128(2),
            denom: T::from_u128(3),
        };
        println!("message: {}", message);
        let ciphertext = crypto_params.encrypt(message.clone());
        println!("ciphertext: {}", ciphertext);
        let decrypted = crypto_params.decrypt(ciphertext);
        assert_eq!(message, decrypted);
    }

    #[test]
    fn public_encrypt_decrypt() {
        let message: Rational<T> = Rational {
            num: T::from_u128(7),
            denom: T::from_u128(3),
        };
        println!("message: {}", message);
        let ciphertext = PUBLIC_PARAMS.encrypt(message.clone());
        println!("ciphertext: {}", ciphertext);
        let decrypted = PUBLIC_PARAMS.decrypt(ciphertext);
        assert_eq!(message, decrypted);
    }

    #[test]
    fn public_encrypt_add_decrypt() {
        let message1: Rational<T> = Rational {
            num: T::from_u128(7),
            denom: T::from_u128(3),
        };

        let message2: Rational<T> = Rational {
            num: T::from_u128(16),
            denom: T::from_u128(5),
        };
        println!("message1: {}", message1);
        println!("message2: {}", message2);
        let ciphertext1 = PUBLIC_PARAMS.encrypt(message1.clone());
        let ciphertext2 = PUBLIC_PARAMS.encrypt(message2.clone());
        let ciphertext_1_plus_2 = ciphertext1 + ciphertext2;
        println!("ciphertext_1_plus_2: {}", ciphertext_1_plus_2);
        let decrypted = PUBLIC_PARAMS.decrypt(ciphertext_1_plus_2);
        assert_eq!(message1 + message2, decrypted);
    }

    #[test]
    fn public_encrypt_mul_decrypt() {
        let message1: Rational<T> = Rational {
            num: T::from_u128(7),
            denom: T::from_u128(3),
        };

        let message2: Rational<T> = Rational {
            num: T::from_u128(16),
            denom: T::from_u128(5),
        };
        println!("message1: {}", message1);
        println!("message2: {}", message2);
        let ciphertext1 = PUBLIC_PARAMS.encrypt(message1.clone());
        let ciphertext2 = PUBLIC_PARAMS.encrypt(message2.clone());
        let ciphertext_1_times_2 = ciphertext1 * ciphertext2;
        println!("ciphertext_1_times_2: {}", ciphertext_1_times_2);
        let decrypted = PUBLIC_PARAMS.decrypt(ciphertext_1_times_2);
        assert_eq!(message1 * message2, decrypted);
    }

    #[test]
    fn public_encrypt_many_add_decrypt() {
        let number_operation = 2;
        let message1: Rational<T> = Rational {
            num: T::from_u128(2),
            denom: T::from_u128(3),
        };

        let message2: Rational<T> = Rational {
            num: T::from_u128(7),
            denom: T::from_u128(4),
        };
        println!("message1: {}", message1);
        println!("message2: {}", message2);
        let ciphertext1 = PUBLIC_PARAMS.encrypt(message1.clone());
        let ciphertext2 = PUBLIC_PARAMS.encrypt(message2.clone());

        let mut clear_result = &message1 + &message2;
        let mut encrypted_result = &ciphertext1 + &ciphertext2;

        for _ in 0..number_operation {
            clear_result = &clear_result + &message1;
            clear_result = &clear_result + &message2;
            encrypted_result = &encrypted_result + &ciphertext1;
            encrypted_result = &encrypted_result + &ciphertext2;
        }

        let decrypted = PUBLIC_PARAMS.decrypt(encrypted_result);
        assert_eq!(clear_result, decrypted);
    }

    // #[test]
    // fn public_encrypt_many_mul_decrypt() {
    // 	let number_operation = 1;
    //     let message1: Rational<T> = Rational {
    //         num: T::from_u128(1),
    //         denom: T::from_u128(1),
    //     };

    // 	let message2: Rational<T> = Rational {
    //         num: T::from_u128(1),
    //         denom: T::from_u128(1),
    //     };
    //     println!("message1: {}", message1);
    // 	println!("message2: {}", message2);
    //     let ciphertext1 = PUBLIC_PARAMS.encrypt(message1.clone());
    // 	let ciphertext2 = PUBLIC_PARAMS.encrypt(message2.clone());

    // 	let mut clear_result = &message1 * &message2;
    // 	let mut encrypted_result = &ciphertext1 * &ciphertext2;

    // 	for _ in 0..number_operation {
    // 	    clear_result = &clear_result *  &message1;
    // 	    clear_result = &clear_result * &message2;
    // 	    encrypted_result = &encrypted_result * &ciphertext1;
    // 	    encrypted_result = &encrypted_result * &ciphertext2;
    // 	}

    //     let decrypted = PUBLIC_PARAMS.decrypt(encrypted_result);
    //     assert_eq!(clear_result, decrypted);
    // }

    //  #[test]
    //     fn public_encrypt_add_decrypt() {
    // 	// let lambda = 256;
    // 	// let d = 2;
    // 	// let eta = d*lambda;
    // 	// let mu = d*d*lambda*9 - eta - 2*lambda - 3; // 9 is log_2(512)
    // 	// let gamma = 2*eta + 3*lambda/2 + mu + 3;
    // 	// assert_eq!(gamma/64, 1);
    // 	// let (p1, p2, p3, p4prime,) = (
    //         //     T::from_u128(7919),
    //         //     T::from_u128(37),
    //         //     T::from_u128(41),
    //         //     T::from_u128(5897),
    //         //     T::from_u128(7759),
    //         // );
    // 	let crypto_param = PublicKeySchemeCryptographicParameters::<T>::new_from_params(256, 2);
    //         let r1 = Rational::<T> {
    //             num: T::from_u128(37),
    //             denom: T::from_u128(7),
    //         };
    //         let r2 = Rational::<T> {
    //             num: T::from_u128(43),
    //             denom: T::from_u128(7),
    //         };
    // 	let r_sum = (&r1 + &r2).reduce();
    //         println!("rational 1: {r1}");
    //         println!("rational 2: {r2}");
    //         let encrypted_1 = crypto_param.encrypt(r1);
    //         let encrypted_2 = crypto_param.encrypt(r2);
    //         println!("encrypted message 1: {encrypted_1}");
    //         println!("encrypted message 2: {encrypted_2}");
    //         let encrypted_sum = encrypted_1 + encrypted_2;
    //         let sum = crypto_param.decrypt(encrypted_sum);
    //         assert_eq!(sum, r_sum);
    // //	assert_eq!((mu, gamma/64), (1, 1));
    //     }

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
