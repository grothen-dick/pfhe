use std::fmt; // used for displaying stuff

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

pub struct Rational {
    num: i64,
    denom: u32,
}

pub struct HenselCode {
    p: u32, // p is assumed to be a prime at this point
    n: i64, // n is in Z/pZ and represents a rational num/denom, where
            // abs(num), abs(denom) <= sqrt((p-1)/2)
}

// pretty print rationals
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

// pretty print hensel codes
impl fmt::Display for HenselCode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} (mod {})", self.n, self.p)
    }
}

pub fn rational_to_hensel_code(p: u32, r: &Rational) -> HenselCode {
    let q = p as i64;
    let n = ((r.num % q) * modular_inverse(r.denom, p)) % q;
    HenselCode { p, n }
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
}
