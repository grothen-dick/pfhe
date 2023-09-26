pub mod bigint;
pub mod crypto_parameters;
pub mod hensel_code;
pub mod rational;
pub mod shared;

use std::{
    clone::Clone,
    fmt, // used for displaying stuff
    ops,
};

#[cfg(test)]
mod tests {
    use super::bigint::BigInt;
    use super::hensel_code::{new_hensel_code, HenselCode};
    use super::rational::Rational;
    use super::shared::DEFAULT_LIMBS;

    const L: usize = DEFAULT_LIMBS;

    #[test]
    fn translates_rational_to_hensel_code() {
        fn simple_tester(r: Rational<L>, p: BigInt) -> () {
            let num = r.num.clone();
            let denom = r.denom.clone();
            let rclone = Rational::<L> {
                num: num.clone(),
                denom: denom.clone(),
            };
            let hc = HenselCode::from((p.clone(), rclone));
            let id_hc = &new_hensel_code(p.clone(), denom).invert();
            let n_hc = new_hensel_code(p.clone(), num);
            assert_eq!(hc.modulus().0, p.0);
            assert_eq!(hc.to_bigint().0, (id_hc * &n_hc).to_bigint().0);
            println!("rational: {} => hensel code: {}", r, hc);
        }

        let p: BigInt = BigInt::<L>::from(37 as u128);

        // positive integer
        let r1 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(1 as u128),
        };
        simple_tester(r1, p.clone());

        // integer inverse
        let r2 = Rational::<L> {
            num: BigInt::<L>::from(1 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(r2, p.clone());

        // general rational
        let r3 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(r3, p);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        fn simple_tester(r: Rational<L>, p: BigInt) -> () {
            let num = r.clone().num;
            let denom = r.clone().denom;
            let hc = HenselCode::from((p.clone(), r.clone()));
            let hcclone = hc.clone();
            let new_r = Rational::<L>::from(hcclone);
            let id_hc = &new_hensel_code(p.clone(), denom).invert();
            let n_hc = new_hensel_code(p.clone(), num);
            assert_eq!(hc.modulus().0, p.0);
            assert_eq!(hc.to_bigint().0, (id_hc * &n_hc).to_bigint().0);
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
        simple_tester(r1, p.clone());

        // integer inverse
        let r2 = Rational::<L> {
            num: BigInt::<L>::from(1 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(r2, p.clone());

        // general rational
        let r3 = Rational::<L> {
            num: BigInt::<L>::from(6 as u128),
            denom: BigInt::<L>::from(8 as u128),
        };
        simple_tester(r3, p);
    }

    #[test]
    fn adds_rationals() {
        fn simple_tester(r1: Rational<L>, r2: Rational<L>) -> () {
            let sum = &r1 + &r2;
            let (n1, n2) = (r1.clone().num, r2.clone().num);
            let (d1, d2) = (r1.clone().denom, r2.clone().denom);
            assert_eq!(sum.denom.0 .0, (d1.clone() * d2.clone()).0 .0);
            assert_eq!(sum.num.0 .0, (n1 * d2 + d1 * n2).0 .0);
            println!("{} + {} = {}", r1, r2, sum);
        }
        // integer addition
        let r11 = Rational::<L> {
            num: BigInt::<L>::from(3 as u128),
            denom: BigInt::<L>::from(1 as u128),
        };
        let r21 = Rational::<L> {
            num: BigInt::<L>::from(1312 as u128),
            denom: BigInt::<L>::from(1 as u128),
        };
        simple_tester(r11, r21);

        // rational addition
        let r12 = Rational::<L> {
            num: BigInt::<L>::from(1337 as u128),
            denom: BigInt::<L>::from(42 as u128),
        };
        let r22 = Rational::<L> {
            num: BigInt::<L>::from(56 as u128),
            denom: BigInt::<L>::from(982 as u128),
        };
        simple_tester(r12, r22);
    }
}
