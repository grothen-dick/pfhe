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

    use super::bigint::BigInt;
    use super::hensel_code::{new_hensel_code, HenselCode};
    use super::rational::Rational;
    use super::shared::DEFAULT_LIMBS;

    const L: usize = DEFAULT_LIMBS;

    #[test]
    fn translates_rational_to_hensel_code() {
        fn simple_tester(r: &Rational<L>, p: &BigInt) -> () {
            let hc = HenselCode::from((p, r));
            let id_hc = &new_hensel_code(p, &(r.denom)).invert();
            let n_hc = new_hensel_code(p, &(r.num));
            assert_eq!(&hc.modulus(), p);
            assert_eq!(hc.to_bigint(), (id_hc * &n_hc).to_bigint());
            println!("rational: {} => hensel code: {}", r, hc);
        }

        let p: BigInt = BigInt::<L>::from(37 as u128);

        // positive integer
        let r1 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(1 as u128),
        };
        simple_tester(&r1, &p);

        // integer inverse
        let r2 = Rational::<L> {
            num: BigInt::<L>::from(1 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(&r2, &p);

        // general rational
        let r3 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(&r3, &p);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        fn simple_tester(r: &Rational<L>, p: &BigInt) -> () {
            let hc = HenselCode::from((p, r));
            let new_r = Rational::<L>::from(&hc);
            let id_hc = new_hensel_code(p, &r.denom).invert();
            let n_hc = new_hensel_code(p, &r.num);
            assert_eq!(&hc.modulus(), p);
            assert_eq!(hc.to_bigint(), (id_hc * n_hc).to_bigint());
            println!(
                "rational: {} => hensel code: {} => rational: {}",
                r, hc, new_r
            );
        }

        let p: BigInt = BigInt::<L>::from(7919 as u128); // thanks wikipedia for this prime

        // positive integer
        let r1 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(1 as u128),
        };
        simple_tester(&r1, &p);

        // integer inverse
        let r2 = Rational::<L> {
            num: BigInt::<L>::from(1 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(&r2, &p);

        // general rational
        let r3 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(&r3, &p);
    }

    #[test]
    fn encrypt_decrypt() {
        let (p1, p2, p3, p4, p5) = (
            BigInt::<L>::from(7919),
            BigInt::<L>::from(37),
            BigInt::<L>::from(41),
            BigInt::<L>::from(5897),
            BigInt::<L>::from(7759),
        );
        let crypto_params: CryptographicParameters<L> =
            CryptographicParameters::<L>::new(p1, p2, p3, p4, p5);
        let message: Rational<L> = Rational {
            num: BigInt::<L>::from(43),
            denom: BigInt::<L>::from(44),
        };
        println!("message: {}", message);
        let ciphertext = crypto_params.encrypt(message.clone());
        println!("ciphertext: {}", ciphertext);
        let decrypted = crypto_params.decrypt(ciphertext);
        assert_eq!(message, decrypted);
    }
}
