// this is for allowing using constants expressions as type parameters
// used mostly for `chinese_remainder`
#![feature(generic_const_exprs)]

pub mod bigint;
pub mod crypto_parameters;
pub mod hensel_code;
pub mod rational;
pub mod shared;

use std::{
    clone::Clone,
    convert::From, // convert hensel code <-> rational
    fmt,           // used for displaying stuff
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
        fn simple_tester(r: Rational, p: BigInt) -> () {
            let num = r.num.clone();
            let denom = r.denom.clone();
            let rclone = Rational { num, denom };
            let hc = HenselCode::from((p, rclone));
            let id_hc = &new_hensel_code(p, denom).invert();
            let n_hc = new_hensel_code(p, num);
            assert_eq!(hc.modulus().0 .0, p.0 .0);
            assert_eq!(hc.to_bigint().0 .0, (id_hc * &n_hc).to_bigint().0 .0);
            println!("rational: {} => hensel code: {}", r, hc);
        }

        let p: BigInt = 37.into();

        // positive integer
        let r1 = Rational {
            num: 6.into(),
            denom: 1.into(),
        };
        simple_tester(r1, p);

        // integer inverse
        let r2 = Rational {
            num: 1.into(),
            denom: 8.into(),
        };
        simple_tester(r2, p);

        // general rational
        let r3 = Rational {
            num: 6.into(),
            denom: 8.into(),
        };
        simple_tester(r3, p);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        fn simple_tester(r: Rational, p: BigInt) -> () {
            let num = r.num;
            let denom = r.denom;
            let rclone = r.clone();
            let hc = HenselCode::from((p, rclone));
            let hcclone = hc.clone();
            let new_r = Rational::from(hcclone);
            let id_hc = &new_hensel_code(p, denom).invert();
            let n_hc = new_hensel_code(p, num);
            assert_eq!(hc.modulus().0 .0, p.0 .0);
            assert_eq!(hc.to_bigint().0 .0, (id_hc * &n_hc).to_bigint().0 .0);
            println!(
                "rational: {} => hensel code: {} => rational: {}",
                r, hc, new_r
            );
        }

        let p: BigInt = 7919.into(); // thanks wikipedia for this prime

        // positive integer
        let r1 = Rational {
            num: 6.into(),
            denom: 1.into(),
        };
        simple_tester(r1, p);

        // integer inverse
        let r2 = Rational {
            num: 1.into(),
            denom: 8.into(),
        };
        simple_tester(r2, p);

        // general rational
        let r3 = Rational {
            num: 6.into(),
            denom: 8.into(),
        };
        simple_tester(r3, p);

        // // negative integer
        // let r4 = Rational { num: -6, denom: 1 };
        // simple_tester(r4, p);

        // // general negative rational
        // let r5 = Rational { num: -7, denom: 8 };
        // simple_tester(r5, p);
    }

    #[test]
    fn adds_rationals() {
        fn simple_tester(r1: Rational, r2: Rational) -> () {
            let sum = &r1 + &r2;
            let (n1, n2) = (r1.num, r2.num);
            let (d1, d2) = (r1.denom, r2.denom);
            assert_eq!(sum.denom.0 .0, (d1 * d2).0 .0);
            assert_eq!(sum.num.0 .0, (n1 * d2 + d1 * n2).0 .0);
            println!("{} + {} = {}", r1, r2, sum);
        }
        // integer addition
        let r11 = Rational {
            num: 3.into(),
            denom: 1.into(),
        };
        let r21 = Rational {
            num: 1312.into(),
            denom: 1.into(),
        };
        simple_tester(r11, r21);

        // rational addition
        let r12 = Rational {
            num: 1337.into(),
            denom: 42.into(),
        };
        let r22 = Rational {
            num: 56.into(),
            denom: 982.into(),
        };
        simple_tester(r12, r22);
    }
}
