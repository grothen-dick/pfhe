use rand::random; // random numbers
use std::fmt; // used for displaying stuff
use std::{
    convert::From,   // convert hensel code <-> rational
    ops::{Add, Mul}, // add and multiply traits
};

pub struct Rational {
    num: i128,
    denom: i128,
}

/// add two `&Rational`s
impl<'a, 'b> Add<&'b Rational> for &'a Rational {
    type Output = Rational;
    fn add(self, other: &'b Rational) -> Rational {
        let (n1, n2) = (self.num, other.num);
        let (d1, d2) = (self.denom, other.denom);
        Rational {
            num: n1 * d2 + n2 * d1,
            denom: d1 * d2,
        }
    }
}

/// add two `Rational`s
impl Add for Rational {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let (n1, n2) = (self.num, other.num);
        let (d1, d2) = (self.denom, other.denom);
        Self {
            num: n1 * d2 + n2 * d1,
            denom: d1 * d2,
        }
    }
}

/// multiply two `&Rational`s
impl<'a, 'b> Mul<&'b Rational> for &'a Rational {
    type Output = Rational;
    fn mul(self, other: &'b Rational) -> Rational {
        let (n1, n2) = (self.num, other.num);
        let (d1, d2) = (self.denom, other.denom);
        Rational {
            num: n1 * n2,
            denom: d1 * d2,
        }
    }
}

/// multiply two `Rational`s
impl Mul for Rational {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let (n1, n2) = (self.num, other.num);
        let (d1, d2) = (self.denom, other.denom);
        Self {
            num: n1 * n2,
            denom: d1 * d2,
        }
    }
}

/// pretty-print Rational
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)   -n_max <= num   <= n_max,
///  ii)  0      <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
impl From<HenselCode> for Rational {
    fn from(hc: HenselCode) -> Self {
        let (p, p_i128, p_f64) = (hc.p, (hc.p as i128), (hc.p as f64));
        let sign_n: i128 = if hc.n < 0 { -1 } else { 1 };
        let abs_n = sign_n * hc.n;

        let n_max = ((p_f64 - 1.0) / 2.0).sqrt() as i128;

        // perform (modified) extended euclidean algorithm on (p, n % p)
        let (mut x0, mut x1): (i128, i128) = (p, i128::try_from(abs_n % p_i128).unwrap());
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

/// add two `&Henselcode`s
impl<'a, 'b> Add<&'b HenselCode> for &'a HenselCode {
    type Output = HenselCode;
    fn add(self, other: &'b HenselCode) -> HenselCode {
        let (n1, n2) = (self.n, other.n);
        let (p1, p2) = (self.p, other.p);

        if p1 != p2 {
            panic!("cannot add '{}' and '{}'", self, other);
        }

        HenselCode { p: p1, n: n1 + n2 }
    }
}

/// add two `Henselcode`s
impl Add for HenselCode {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let (n1, n2) = (self.n, other.n);
        let (p1, p2) = (self.p, other.p);

        if p1 != p2 {
            panic!("cannot add '{}' and '{}'", self, other);
        }

        Self { p: p1, n: n1 + n2 }
    }
}

/// multiply two `&Henselcode`s
impl<'a, 'b> Mul<&'b HenselCode> for &'a HenselCode {
    type Output = HenselCode;
    fn mul(self, other: &'b HenselCode) -> HenselCode {
        let (n1, n2) = (self.n, other.n);
        let (p1, p2) = (self.p, other.p);

        if p1 != p2 {
            panic!("cannot multiply '{}' and '{}'", self, other);
        }
        HenselCode { p: p1, n: n1 * n2 }
    }
}

/// multiply two `Henselcode`s
impl Mul for HenselCode {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let (n1, n2) = (self.n, other.n);
        let (p1, p2) = (self.p, other.p);

        if p1 != p2 {
            panic!("cannot multiply '{}' and '{}'", self, other);
        }
        Self { p: p1, n: n1 * n2 }
    }
}

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
        let n = ((r.num % p) * id) % p;
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

    /// given `n_{1}, n_{2}, n_{3}`, compute N such that `N = n_{i} (mod p_{i})`
    fn chinese_remainder(&self, n1: i128, n2: i128, n3: i128) -> HenselCode {
        let (p1, p2, p3) = (self.p1, self.p2, self.p3);
        // we want to solve:
        // ap1p2 + bp1p3 + cp2p3 = n1 (mod p1) = n2 (mod p2) = n3 (mod p3)
        // so we have:
        //     cp2p3 = n1 (mod p1)
        //     bp1p3 = n2 (mod p2)
        //     ap1p2 = n3 (mod p3)
        let (q12, q13, q23) = (p1 * p2, p1 * p3, p2 * p3);
        let (i12, _) = modular_inverses(q12, p3);
        let (i13, _) = modular_inverses(q13, p2);
        let (i23, _) = modular_inverses(q23, p1);
        let n = (n3 * i12 * q12 + n2 * i13 * q13 + n1 * i23 * q23) % (p1 * p2 * p3);
        HenselCode { p: p1 * p2 * p3, n }
    }
    pub fn encode(&self, m: i128) -> HenselCode {
        // TODO: use correct bounds for the variables
        let g = self.public_key();
        let delta_max = g / self.p4;
        let s2 = random::<i128>() % self.p2;
        let s3 = random::<i128>() % self.p3;
        let delta = random::<i128>() % delta_max;

        let rm = Rational { num: m, denom: 1 };
        let s1 = Rational {
            num: random::<i128>(),
            denom: 1,
        };
        let dp4 = HenselCode {
            p: g,
            n: delta * self.p4,
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
            assert_eq!(hc.n, (num * id) % p);
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
            assert_eq!(hc.n, (num * id) % p);
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
}
