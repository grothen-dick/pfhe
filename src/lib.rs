use rand::random; // random numbers
use std::fmt; // used for displaying stuff
use std::ops::{Add, Mul}; // add and multiply traits

pub struct Rational {
    num: i128,
    denom: i128,
}

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

pub struct HenselCode {
    p: i128, // p is assumed to be a prime at this point
    n: i128, // n is in Z/pZ and represents a rational num/denom, where
             // abs(num) <= sqrt((p-1)/2), denom <= 2*sqrt((p-1)/2)
}

/// pretty-print HenselCode
impl fmt::Display for HenselCode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} (mod {})", self.n, self.p)
    }
}

/// compute `y` such that y*m = 1 (mod n)
/// (assuming that there exists such a `y`)
pub fn modular_inverses(n: i128, m: i128) -> (i128, i128) {
    // TODO: return a pair (i_n: i128, i_m: i128) such that
    // i_n * n + i_m * m = 1
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

/// given a prime `p` and a rational `r = num/denom`, where p does not divide denom,
/// returns `r (mod p)`
pub fn rational_to_hensel_code(p: i128, r: &Rational) -> HenselCode {
    let (i_d, _) = modular_inverses(r.denom, p); // REVIEW
    let n = ((r.num % p) * i_d) % p;
    HenselCode { p, n }
}

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)   -n_max <= num   <= n_max,
///  ii)  0      <= denom <= 2*n_max,
///  iii) hc = num/denom (mod p)
pub fn hensel_code_to_rational(hc: &HenselCode) -> Rational {
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
    fn chinese_remainder(self, n1: i128, n2: i128, n3: i128) -> i128 {
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
        (n3 * i12 * q12 + n2 * i13 * q13 + n1 * i23 * q23) % (p1 * p2 * p3)
    }
    pub fn encode(&self, m: i128) -> i128 {
        // TODO: use correct bounds for the variables
        let p4 = self.p4;
        let g = self.public_key();
        let delta_max = g / p4;
        let s1 = random::<i128>();
        let s2 = random::<i128>() % i128::from(self.p2);
        let s3 = random::<i128>() % i128::from(self.p3);
        let delta = random::<i128>() % delta_max;
        // (s1 * () + delta * p4) % g
        return 0;
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
            let num = r.num;
            let denom = r.denom;
            let hc = rational_to_hensel_code(p, &r);
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
            let hc = rational_to_hensel_code(p, &r);
            let new_r = hensel_code_to_rational(&hc);
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
}
