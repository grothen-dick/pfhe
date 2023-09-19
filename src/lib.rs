#[macro_use]
extern crate impl_ops;

use std::fmt; // used for displaying stuff
use std::{
    convert::From, // convert hensel code <-> rational
    ops,           // add and multiply traits
};

pub struct Rational {
    num: i128,
    denom: i128,
}

impl_op_ex!(+ |r1: &Rational, r2: &Rational| -> Rational {
Rational {
    num: r1.num * r2.denom + r2.num * r1.denom,
    denom: r1.denom * r2.denom,
    }
});

impl_op_ex!(*|r1: &Rational, r2: &Rational| -> Rational {
    Rational {
        num: r1.num * r2.num,
        denom: r1.denom * r2.denom,
    }
});

/// pretty-print Rational
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)   0 <= num   <= n_max,
///  ii)  0 <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
impl From<HenselCode> for Rational {
    fn from(hc: HenselCode) -> Self {
        let (p, p_i128, p_f64) = (hc.p, (hc.p as i128), (hc.p as f64));
        let sign_n: i128 = if hc.n < 0 { -1 } else { 1 };
        let abs_n = sign_n * hc.n;

        let n_max = ((p_f64 - 1.0) / 2.0).sqrt() as i128;

        // perform (modified) extended euclidean algorithm on (p, n % p)
        let (mut x0, mut x1): (i128, i128) = (p, abs_n.rem_euclid(p_i128));
        let (mut y0, mut y1): (i128, i128) = (0, 1);
        while x0 > n_max {
            let q = x0 / x1;
            (x0, x1) = (x1, x0 - q * x1);
            (y0, y1) = (y1, y0 - i128::from(q) * y1);
        }
        let (mut num, mut denom): (i128, i128) = (sign_n * (x0 as i128), y0);

        // make sure `denom` is positive
        if denom < 0 {
            (num, denom) = (-num, -denom);
        }
        // return `num/denom`
        Rational {
            num,
            denom: i128::try_from(denom).unwrap(),
        }
    }
}

pub struct HenselCode {
    p: i128, // p is assumed to be a prime at this point
    n: i128, // n is in Z/pZ and represents a rational num/denom, where
             // abs(num) <= sqrt((p-1)/2), denom <= 2*sqrt((p-1)/2)
}

impl HenselCode {
    pub fn chinese_remainder(hc1: Self, hc2: Self) -> Self {
        let (p1, n1) = (hc1.p, hc1.n);
        let (p2, n2) = (hc2.p, hc2.n);
        let (i1, i2) = modular_inverses(p1, p2);
        HenselCode {
            p: p1 * p2,
            n: (p1 * i1 * n2 + p2 * i2 * n1).rem_euclid(p1 * p2),
        }
    }
}

impl_op_ex!(+ |hc1: &HenselCode, hc2: &HenselCode| -> HenselCode {
    if hc1.p != hc2.p {
            panic!("cannot add '{}' and '{}'", hc1, hc2);
    }
    HenselCode {p: hc1.p, n: hc1.n + hc2.n}
});

impl_op_ex!(*|hc1: &HenselCode, hc2: &HenselCode| -> HenselCode {
    if hc1.p != hc2.p {
        panic!("cannot multiply '{}' and '{}'", hc1, hc2);
    }
    HenselCode {
        p: hc1.p,
        n: hc1.n * hc2.n,
    }
});

/// pretty-print HenselCode
impl fmt::Display for HenselCode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} (mod {})", self.n, self.p)
    }
}

/// given a prime `p` and a rational `r = num/denom`, where p does not divide denom,
/// returns `r (mod p)`
impl From<(i128, Rational)> for HenselCode {
    fn from(params: (i128, Rational)) -> Self {
        let (p, r) = params;
        let (id, _) = modular_inverses(r.denom, p);
        let n = ((r.num.rem_euclid(p)) * id).rem_euclid(p);
        HenselCode { p, n }
    }
}

/// return a pair (i_n, i_m) such that
/// i_n * n + i_m * m = 1
pub fn modular_inverses(n: i128, m: i128) -> (i128, i128) {
    if n < m {
        let (i_m, i_n) = modular_inverses(m, n);
        return (i_n, i_m);
    }
    let (mut x0, mut x1) = (n, m);
    let (mut y0, mut y1): (i128, i128) = (1, 0);
    let (mut z0, mut z1): (i128, i128) = (0, 1);
    while x1 > 0 {
        let q = x0 / x1;
        (x0, x1) = (x1, x0 - q * x1);
        (y0, y1) = (y1, y0 - q * y1);
        (z0, z1) = (z1, z0 - q * z1);
    }
    (y0, z0)
}

pub struct CryptographicParameters {
    p1: i128,
    p2: i128,
    p3: i128,
    p4: i128,
    p5: i128,
}

impl CryptographicParameters {
    /// return the product of the 5 primes used as crypto parameters
    pub fn public_key(&self) -> i128 {
        self.p1 * self.p2 * self.p3 * self.p4 * self.p5
    }

    fn chinese_remainder(&self, n1: i128, n2: i128, n3: i128) -> HenselCode {
        let hc1 = HenselCode { p: self.p1, n: n1 };
        let hc2 = HenselCode { p: self.p2, n: n2 };
        let hc3 = HenselCode { p: self.p3, n: n3 };
        HenselCode::chinese_remainder(HenselCode::chinese_remainder(hc1, hc2), hc3)
    }

    pub fn encode(&self, m: i128) -> HenselCode {
        // TODO: use correct bounds for the variables
        let g = self.public_key();
        // let delta_max = g / self.p4;
        // let s2 = random::<i128>() % self.p2;
        // let s3 = random::<i128>() % self.p3;
        // let delta = random::<i128>() % delta_max;
        let s2 = 0; // FIXME
        let s3 = 0; // FIXME
        let delta = 0; // FIXME

        let dp4 = HenselCode {
            p: g,
            n: delta * self.p4,
        };

        let rm = Rational { num: m, denom: 1 };
        let s1 = Rational {
            // num: random::<i128>(),
            num: 0, // FIXME
            denom: 1,
        };
        // intermediary step with a Farey fraction
        let rational_term = s1 * Rational::from(self.chinese_remainder(0, s2, s3)) + rm;
        // returns a HenselCode
        return HenselCode::from((g, rational_term)) + dp4;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_modular_inverses() {
        let inv = modular_inverses(37, 5);
        assert_eq!(inv, (-2, 15)); // -2*37 + 15*5 = 1
    }

    #[test]
    fn finds_modular_inverses_big() {
        let inv = modular_inverses(374533, 52343);
        assert_eq!(inv, (-19870, 142177));
    }

    #[test]
    fn translates_rational_to_hensel_code() {
        fn simple_tester(r: Rational, p: i128) -> () {
            let num = r.num.clone();
            let denom = r.denom.clone();
            let rclone = Rational { num, denom };
            let hc = HenselCode::from((p, rclone));
            let (id, _) = modular_inverses(denom, p);
            assert_eq!(hc.p, p);
            assert_eq!(hc.n, (num * id).rem_euclid(p));
            println!("rational: {} => hensel code: {}", r, hc);
        }

        let p: i128 = 37;

        // positive integer
        let r1 = Rational { num: 6, denom: 1 };
        simple_tester(r1, p);

        // integer inverse
        let r2 = Rational { num: 1, denom: 8 };
        simple_tester(r2, p);

        // general rational
        let r3 = Rational { num: 6, denom: 8 };
        simple_tester(r3, p);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        fn simple_tester(r: Rational, p: i128) -> () {
            let num = r.num;
            let denom = r.denom;
            let rclone = Rational { num, denom };
            let hc = HenselCode::from((p, rclone));
            let hcclone = HenselCode { p, n: hc.n };
            let new_r = Rational::from(hcclone);
            let (id, _) = modular_inverses(denom, p);
            assert_eq!(hc.p, p);
            assert_eq!(hc.n, (num * id).rem_euclid(p));
            println!(
                "rational: {} => hensel code: {} => rational: {}",
                r, hc, new_r
            );
        }

        let p: i128 = 7919; // thanks wikipedia for this prime

        // positive integer
        let r1 = Rational { num: 6, denom: 1 };
        simple_tester(r1, p);

        // integer inverse
        let r2 = Rational { num: 1, denom: 8 };
        simple_tester(r2, p);

        // general rational
        let r3 = Rational { num: 6, denom: 8 };
        simple_tester(r3, p);

        // negative integer
        let r4 = Rational { num: -6, denom: 1 };
        simple_tester(r4, p);

        // general negative rational
        let r5 = Rational { num: -7, denom: 8 };
        simple_tester(r5, p);
    }

    #[test]
    fn adds_rationals() {
        fn simple_tester(r1: Rational, r2: Rational) -> () {
            let sum = &r1 + &r2;
            let (n1, n2) = (r1.num, r2.num);
            let (d1, d2) = (r1.denom, r2.denom);
            assert_eq!(sum.denom, d1 * d2);
            assert_eq!(sum.num, n1 * d2 + n2 * d1);
            println!("{} + {} = {}", r1, r2, sum);
        }
        // integer addition
        let r11 = Rational { num: 3, denom: 1 };
        let r21 = Rational {
            num: 1312,
            denom: 1,
        };
        simple_tester(r11, r21);

        // rational addition
        let r12 = Rational {
            num: -1337,
            denom: 42,
        };
        let r22 = Rational {
            num: 56,
            denom: 982,
        };
        simple_tester(r12, r22);
    }

    #[test]
    fn chinese_remainder() {
        let (p1, p2, p3, p4, p5) = (4919, 7, 11, 13, 17);
        let crypto_param = CryptographicParameters { p1, p2, p3, p4, p5 };
        let (n1, n2, n3) = (38, 2, 1);
        let result = crypto_param.chinese_remainder(n1, n2, n3);
        // for 0 reason the usual 'modulo' operator `%` has a range of [-p, p]
        // for `(...) % p` instead of the normal range [0, p]. The sensible version
        // is instead written `(...).rem_euclid(p)`.
        assert_eq!(result.n.rem_euclid(4919), n1);
        assert_eq!(result.n.rem_euclid(7), n2);
        assert_eq!(result.n.rem_euclid(11), n3);
        let hc1 = HenselCode { p: p1, n: n1 };
        let hc2 = HenselCode { p: p2, n: n2 };
        let hc3 = HenselCode { p: p3, n: n3 };
        let hc12 = HenselCode::chinese_remainder(hc1, hc2);
        let hc = HenselCode::chinese_remainder(hc12, hc3);
        assert_eq!(result.n.rem_euclid(p1 * p2 * p3), hc.n);
        println!("{} : {}", hc, result);
    }
}
