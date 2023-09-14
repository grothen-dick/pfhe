use std::fmt; // used for displaying stuff

pub struct Rational {
    num: i64,
    denom: u32,
}

pub struct HenselCode {
    p: u32, // p is assumed to be a prime at this point
    n: i64, // n is in Z/pZ and represents a rational num/denom, where
            // abs(num), abs(denom) <= sqrt((p-1)/2)
}

/// pretty print Rational
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

/// pretty print HenselCode
impl fmt::Display for HenselCode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} (mod {})", self.n, self.p)
    }
}

/// compute `y` such that y*num2 = 1 (mod num1)
/// (assuming that there exists such a `y`)
pub fn modular_inverse(num1: u32, num2: u32) -> i64 {
    let (mut x0, mut x1) = if num1 > num2 {
        (num1, num2)
    } else {
        (num2, num1)
    };
    let (mut y0, mut y1): (i64, i64) = (0, 1);
    while x1 > 0 {
        let q = x0 / x1;
        (x0, x1) = (x1, x0 - q * x1);
        (y0, y1) = (y1, y0 - i64::from(q) * y1);
    }
    y0
}

/// given a prime `p` and a rational `r = num/denom`, where p does not divide denom,
/// returns `r (mod p)`
pub fn rational_to_hensel_code(p: u32, r: &Rational) -> HenselCode {
    let q = p as i64;
    let n = ((r.num % q) * modular_inverse(r.denom, p)) % q;
    HenselCode { p, n }
}

/// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
/// num/denom where:
///  i)  abs(num), abs(denom) <= n_max,
///  ii) hc = num/denom (mod p)
pub fn hensel_code_to_rational(hc: &HenselCode) -> Rational {
    let (p, p_i64, p_f64) = (hc.p, (hc.p as i64), (hc.p as f64));
    let n_max = ((p_f64 - 1.0) / 2.0).sqrt() as u32;
    // perform (modified) extended euclidean algorithm on (p, hc.n % p)
    let (mut x0, mut x1): (u32, u32) = (p, u32::try_from(hc.n % p_i64).unwrap());
    let (mut y0, mut y1): (i64, i64) = (0, 1);
    while x0 > n_max {
        let q = x0 / x1;
        (x0, x1) = (x1, x0 - q * x1);
        (y0, y1) = (y1, y0 - i64::from(q) * y1);
    }
    let (mut num, mut denom): (i64, i64) = (x0 as i64, y0 as i64);

    // make sure `denom` is positive
    if denom < 0 {
        (num, denom) = (-num, -denom);
    }
    // return `num/denom`
    Rational {
        num,
        denom: u32::try_from(denom).unwrap(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_finds_modular_inverse() {
        let inv_x = modular_inverse(37, 5);
        assert_eq!(inv_x, 15);
    }

    #[test]
    fn it_finds_modular_inverse_wrong_order() {
        let inv_x = modular_inverse(374533, 52343);
        assert_eq!(inv_x, 142177);
    }

    #[test]
    fn translates_rational_to_hensel_code() {
        let p: u32 = 37;
        let q = p as i64;

        let r1 = Rational { num: 6, denom: 1 };
        let hc1 = rational_to_hensel_code(p, &r1);
        assert_eq!(hc1.p, p);
        assert_eq!(hc1.n, r1.num % q);
        println!("rational: {} => hensel code: {}", r1, hc1);

        let r2 = Rational { num: 1, denom: 8 };
        let hc2 = rational_to_hensel_code(p, &r2);
        assert_eq!(hc2.p, p);
        assert_eq!(hc2.n, modular_inverse(p, r2.denom));
        println!("rational: {} => hensel code: {}", r2, hc2);

        let r3 = Rational { num: 6, denom: 8 };
        let hc3 = rational_to_hensel_code(p, &r3);
        assert_eq!(hc3.p, p);
        assert_eq!(hc3.n, (r3.num * modular_inverse(p, r3.denom)) % q);
        println!("rational: {} => hensel code: {}", r3, hc3);
    }

    #[test]
    fn translates_rational_to_hc_and_back() {
        let p: u32 = 7919; // thanks wikipedia for this prime
        let q = p as i64;

        let r1 = Rational { num: 6, denom: 1 };
        let hc1 = rational_to_hensel_code(p, &r1);
        let new_r1 = hensel_code_to_rational(&hc1);
        assert_eq!(hc1.p, p);
        assert_eq!(hc1.n, r1.num % q);
        println!(
            "rational: {} => hensel code: {} => rational: {}",
            r1, hc1, new_r1
        );

        let r2 = Rational { num: 1, denom: 8 };
        let hc2 = rational_to_hensel_code(p, &r2);
        let new_r2 = hensel_code_to_rational(&hc2);
        assert_eq!(hc2.p, p);
        assert_eq!(hc2.n, modular_inverse(p, r2.denom));
        println!(
            "rational: {} => hensel code: {} => rational: {}",
            r2, hc2, new_r2
        );

        let r3 = Rational { num: 6, denom: 8 };
        let hc3 = rational_to_hensel_code(p, &r3);
        let new_r3 = hensel_code_to_rational(&hc3);
        assert_eq!(hc3.p, p);
        assert_eq!(hc3.n, (r3.num * modular_inverse(p, r3.denom)) % q);
        println!(
            "rational: {} => hensel code: {} => rational: {}",
            r3, hc3, new_r3
        );
    }
}
