// this is for allowing using constants expressions as type parameters
// used mostly for `chinese_remainder`
#![feature(generic_const_exprs)]

mod crypto_parameters;

use std::{
    clone::Clone,
    convert::From, // convert hensel code <-> rational
    fmt,           // used for displaying stuff
    ops,
};

use crypto_bigint::U256;
/// default size of a BigInt (in LIMBS)
const DEFAULT_LIMBS: usize = U256::LIMBS;

/// A trait to retrieve the current size of a const-generic struct
trait Bounded {
    const L: usize;

    fn size(&self) -> usize {
        Self::L
    }
}

mod big_int {
    use super::{
        fmt,
        ops::{Add, Div, Mul, Rem, Sub},
        Bounded, Clone, From, DEFAULT_LIMBS,
    };
    use crypto_bigint::{rand_core::OsRng, NonZero, RandomMod, Uint, Wrapping, Zero, U128};

    /// Simple wrapper to abstract details away from crypto_bigint library.
    /// We simply want to be able to:
    /// - do arithmetics with BigInt
    /// - display it with println! for debug purposes
    /// - create a BigInt from a regular integer
    /// - generate a random BigInt
    /// - change the size of BigInt
    #[derive(PartialEq, PartialOrd)]
    pub struct BigInt<const L: usize = DEFAULT_LIMBS>(pub Wrapping<Uint<L>>);

    impl<const L: usize> Bounded for BigInt<L> {
        const L: usize = L;
    }

    impl<const L: usize> BigInt<L> {
        /// create a random BigInt modulo `modulus`
        pub fn random_mod(modulus: &BigInt<L>) -> BigInt<L> {
            BigInt(Wrapping(Uint::<L>::random_mod(
                &mut OsRng,
                &NonZero::new(modulus.0 .0).unwrap(),
            )))
        }

        /// wrap a Uint into a BigInt
        pub fn new(n: Uint<L>) -> BigInt<L> {
            BigInt(Wrapping(n))
        }

        pub fn resize<const Lnew: usize>(&self) -> BigInt<Lnew> {
            BigInt::new(self.0 .0.resize::<Lnew>())
        }

        /// compute square root of BigInt
        pub fn sqrt(self) -> BigInt<L> {
            BigInt(Wrapping(self.0 .0.sqrt_vartime()))
        }

        /// compute gcd of two BigInt, the good-old Euclid way
        pub fn gcd(b1: &BigInt<L>, b2: &BigInt<L>) -> BigInt<L> {
            if b1 < b2 {
                return Self::gcd(b2, b1);
            }
            let (mut x0, mut x1) = (b1.clone(), b2.clone());
            while x1.0 .0.is_zero().unwrap_u8() == 0 {
                (x1, x0) = (&x0 % &x1, x1);
            }
            return x0;
        }
    }
    // clone a BgInt
    impl<const L: usize> Clone for BigInt<L> {
        fn clone(&self) -> Self {
            BigInt::new(self.0 .0.clone())
        }
    }

    // creates a BigInt from a regular integer
    impl<const L: usize> From<u128> for BigInt<L> {
        fn from(k: u128) -> BigInt<L> {
            BigInt(Wrapping(Uint::<L>::from(k)))
        }
    }

    impl<const L: usize> From<&u128> for BigInt<L> {
        fn from(k: &u128) -> BigInt<L> {
            BigInt(Wrapping(Uint::<L>::from(k.clone())))
        }
    }

    // display a BigInt
    impl<const L: usize> fmt::Display for BigInt<L> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{}", self.0 .0.to_string())
        }
    }

    // implement add, sub, mul, div, rem for &BigInt
    impl<'a, 'b, const L: usize> Add<&'b BigInt<L>> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn add(self, other: &'b BigInt<L>) -> BigInt<L> {
            let (b1, b2) = (self.0, other.0);
            BigInt(b1 + b2)
        }
    }
    impl<'a, 'b, const L: usize> Sub<&'b BigInt<L>> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn sub(self, other: &'b BigInt<L>) -> BigInt<L> {
            let (b1, b2) = (self.0, other.0);
            BigInt(b1 - b2)
        }
    }
    impl<'a, 'b, const L1: usize, const L2: usize> Mul<&'b BigInt<L2>> for &'a BigInt<L1>
    where
        BigInt<{ L1 + L2 }>: Sized,
    {
        type Output = BigInt<{ L1 + L2 }>;
        fn mul(self, other: &'b BigInt<L2>) -> BigInt<{ L1 + L2 }> {
            let (b1, b2) = (
                self.resize::<{ L1 + L2 }>().0,
                other.resize::<{ L1 + L2 }>().0,
            );
            BigInt::<{ L1 + L2 }>(b1 * b2)
        }
    }
    impl<'a, 'b, const L: usize> Div<&'b BigInt<L>> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn div(self, other: &'b BigInt<L>) -> BigInt<L> {
            let (b1, b2) = (self.0, NonZero::new(other.0 .0).unwrap());
            BigInt(b1 / b2)
        }
    }
    impl<'a, 'b, const L: usize> Rem<&'b BigInt<L>> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn rem(self, other: &'b BigInt<L>) -> BigInt<L> {
            let (b1, b2) = (self.0, NonZero::new(other.0 .0).unwrap());
            BigInt(b1 % b2)
        }
    }

    // implement operations for BigInt
    impl<const L: usize> Add<BigInt<L>> for BigInt<L> {
        type Output = BigInt<L>;
        fn add(self, other: BigInt<L>) -> BigInt<L> {
            &self + &other
        }
    }
    impl<const L: usize> Sub<BigInt<L>> for BigInt<L> {
        type Output = BigInt<L>;
        fn sub(self, other: BigInt<L>) -> BigInt<L> {
            &self - &other
        }
    }
    impl<const L1: usize, const L2: usize> Mul<BigInt<L2>> for BigInt<L1>
    where
        BigInt<{ L1 + L2 }>: Sized,
    {
        type Output = BigInt<{ L1 + L2 }>;
        fn mul(self, other: BigInt<L2>) -> BigInt<{ L1 + L2 }> {
            &self * &other
        }
    }
    impl<const L: usize> Div<BigInt<L>> for BigInt<L> {
        type Output = BigInt<L>;
        fn div(self, other: BigInt<L>) -> BigInt<L> {
            &self / &other
        }
    }
    impl<const L: usize> Rem<BigInt<L>> for BigInt<L> {
        type Output = BigInt<L>;
        fn rem(self, other: BigInt<L>) -> BigInt<L> {
            &self % &other
        }
    }

    // implement arithmetic operations with u128
    impl<'a, const L: usize> Add<u128> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn add(self, other: u128) -> BigInt<L> {
            self + &BigInt::from(other)
        }
    }
    impl<'a, const L: usize> Sub<u128> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn sub(self, other: u128) -> BigInt<L> {
            self - &BigInt::from(other)
        }
    }
    impl<'a, const L: usize> Mul<u128> for &'a BigInt<L>
    where
        BigInt<{ L + U128::LIMBS }>: Sized,
        BigInt<{ U128::LIMBS }>: Sized,
    {
        type Output = BigInt<{ L + U128::LIMBS }>;
        fn mul(self, other: u128) -> BigInt<{ L + U128::LIMBS }> {
            self * &BigInt::<{ U128::LIMBS }>::from(other)
        }
    }
    impl<'a, const L: usize> Div<u128> for &'a BigInt<L> {
        type Output = BigInt<L>;
        fn div(self, other: u128) -> BigInt<L> {
            self / &BigInt::from(other)
        }
    }
}

mod rational {

    use super::{
        big_int::BigInt,
        fmt,
        hensel_code::HenselCode,
        ops::{Add, Mul},
        Bounded, Clone, DEFAULT_LIMBS,
    };

    pub struct Rational<const L: usize = DEFAULT_LIMBS> {
        pub num: BigInt<L>,
        pub denom: BigInt<L>,
    }

    impl<const L: usize> Bounded for Rational<L> {
        const L: usize = L;
    }

    impl<const L: usize> Clone for Rational<L> {
        fn clone(&self) -> Self {
            Rational {
                num: self.num.clone(),
                denom: self.denom.clone(),
            }
        }
    }

    /// we want to simplify a rational
    impl<const L: usize> Rational<L> {
        pub fn reduce(&self) -> Rational<L> {
            let gcd = BigInt::gcd(&self.num, &self.denom);
            let num = &self.num / &gcd;
            let denom = &self.denom / &gcd;
            Rational { num, denom }
        }
        pub fn resize<const Lnew: usize>(&self) -> Rational<Lnew> {
            Rational {
                num: self.num.resize::<Lnew>(),
                denom: self.denom.resize::<Lnew>(),
            }
        }
    }

    /// Add two &Rational
    impl<'a, 'b, const L: usize> Add<&'b Rational<L>> for &'a Rational<L>
    where
        [(); L + L]:,
        BigInt<{ L + L }>: Sized,
        Rational<{ L + L }>: Sized,
    {
        type Output = Rational<{ L + L }>;
        fn add(self, other: &'b Rational<L>) -> Rational<{ L + L }> {
            let (r1, r2) = (self, other);
            let num1 = &r1.num * &r2.denom;
            let num2 = &r1.denom * &r2.num; // mystical shit: for Rust, L + L != L + L
            let num = num1 + num2;
            let denom: BigInt<{ L + L }> = &r1.denom * &r2.denom;
            Rational { num, denom }.reduce()
        }
    }

    /// Multiply two &Rational
    impl<'a, 'b, const L1: usize, const L2: usize> Mul<&'b Rational<L2>> for &'a Rational<L1>
    where
        Rational<{ L1 + L2 }>: Sized,
    {
        type Output = Rational<{ L1 + L2 }>;
        fn mul(self, other: &'b Rational<L2>) -> Rational<{ L1 + L2 }> {
            let (r1, r2) = (self, other);
            Rational {
                num: &r1.num * &r2.num,
                denom: &r1.denom * &r2.denom,
            }
            .reduce()
        }
    }

    /// Add two Rational
    impl<const L: usize> Add<Rational<L>> for Rational<L>
    where
        [(); L + L]:,
        Rational<{ L + L }>: Sized,
    {
        type Output = Rational<{ L + L }>;
        fn add(self, other: Rational<L>) -> Rational<{ L + L }> {
            &self + &other
        }
    }

    /// Multiply two Rational
    impl<const L1: usize, const L2: usize> Mul<Rational<L2>> for Rational<L1>
    where
        [(); L1 + L2]:,
        Rational<{ L1 + L2 }>: Sized,
    {
        type Output = Rational<{ L1 + L2 }>;
        fn mul(self, other: Rational<L2>) -> Rational<{ L1 + L2 }> {
            let (r1, r2) = (self, other);
            Rational {
                num: &r1.num * &r2.num,
                denom: &r1.denom * &r2.denom,
            }
            .reduce()
        }
    }

    /// pretty-print Rational
    impl<const L: usize> fmt::Display for Rational<L> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{}/{}", self.num, self.denom)
        }
    }

    /// given an element `hc` of Z/pZ, compute n_max = floor(sqrt((p-1)/2)) and return a rational
    /// num/denom where:
    ///  i)   0 <= num   <= n_max,
    ///  ii)  0 <= denom <= 2*n_max,
    ///  iii) hc = num/denom (mod p)
    impl<const L: usize> From<HenselCode<L>> for Rational<L>
    where
        [(); L]:,
        [(); L + L]:,
        BigInt<{ L + L }>: Sized,
    {
        fn from(hc: HenselCode<L>) -> Self {
            let n_max = (&(&hc.modulus() - 1) / 2).sqrt();

            // perform (modified) extended euclidean algorithm on (g, n % g)
            let (mut x0, mut x1) = (hc.modulus(), hc.to_big_int());
            let (mut y0, mut y1) = (BigInt::<L>::from(0), BigInt::<L>::from(1));
            while x0 > n_max {
                let q = (&x0 / &x1).resize::<L>();
                let new_x0 = x1.clone();
                let big_x1 = x0.resize::<{ L + L }>() - (&q * &x1);
                let new_y0 = y1.clone();
                let big_y1 = y0.resize::<{ L + L }>() - (&q * &y1);
                (x0, x1) = (new_x0, big_x1.resize::<L>());
                (y0, y1) = (new_y0, big_y1.resize::<L>());
            }

            Rational { num: x0, denom: y0 }
        }
    }
}

mod hensel_code {
    use super::{
        big_int::BigInt,
        fmt,
        ops::{Add, Mul},
        rational::Rational,
        Bounded, Clone, DEFAULT_LIMBS,
    };
    use crypto_bigint::modular::runtime_mod::{DynResidue, DynResidueParams};

    // the operation `chinese_remainder` changes the size of the modulus, so we need to track it using a const generics LIMBS
    pub struct HenselCode<const L: usize = DEFAULT_LIMBS> {
        params: DynResidueParams<L>,
        res: DynResidue<L>, // internal variable that stores the residue
    }
    impl<const L: usize> HenselCode<L> {
        /// return a BigInt n, with residue `res` mod `params.modulus()`
        pub fn to_big_int(&self) -> BigInt<L> {
            BigInt::new(self.res.retrieve())
        }

        /// return the modulus
        pub fn modulus(&self) -> BigInt<L> {
            BigInt::new(*self.params.modulus())
        }
    }

    impl<const L: usize> Bounded for HenselCode<L> {
        const L: usize = L;
    }

    impl<const L: usize> Clone for HenselCode<L> {
        fn clone(&self) -> Self {
            HenselCode {
                params: self.params.clone(),
                res: self.res.clone(),
            }
        }
    }

    /// Create an HenselCode from two BigInt
    pub fn new_hensel_code<const Lg: usize, const Ln: usize>(
        g: BigInt<Lg>,
        n: BigInt<Ln>,
    ) -> HenselCode<Lg> {
        let params = DynResidueParams::new(&g.0 .0);
        let res = DynResidue::new(&(&n.resize::<Lg>() % &g).0 .0, params);
        HenselCode { params, res }
    }

    impl<const L: usize> HenselCode<L> {
        pub fn invert(&self) -> HenselCode<L> {
            HenselCode {
                params: self.params,
                res: self.res.invert().0,
            }
        }
        pub fn chinese_remainder<const L1: usize, const L2: usize>(
            hc1: HenselCode<L1>,
            hc2: HenselCode<L2>,
        ) -> HenselCode<{ L1 + L2 }> {
            let (g1, n1) = (
                hc1.modulus().resize::<{ L1 + L2 }>(),
                hc1.to_big_int().resize::<{ L1 + L2 }>(),
            );
            let (g2, n2) = (
                hc2.modulus().resize::<{ L1 + L2 }>(),
                hc2.to_big_int().resize::<{ L1 + L2 }>(),
            );
            let g12 = &hc1.modulus() * &hc2.modulus();
            let (residue_params1, residue_params2, residue_params) = (
                DynResidueParams::new(&g1.0 .0),
                DynResidueParams::new(&g2.0 .0),
                DynResidueParams::new(&g12.0 .0),
            );
            let (mut res_g1, mut res_g2) = (
                DynResidue::new(&g1.0 .0, residue_params2),
                DynResidue::new(&g2.0 .0, residue_params1),
            );
            // i1*g1 = 1 (mod g2), i2*g2 = 1 (mod g1)
            // we need to convert i1 and i2 to a residue mod g1*g2
            let (i1, i2) = (
                DynResidue::new(&res_g1.invert().0.retrieve(), residue_params.clone()),
                DynResidue::new(&res_g2.invert().0.retrieve(), residue_params.clone()),
            );

            // change modulus g1 -> g1*g2 and g2 -> g1*g2
            (res_g1, res_g2) = (
                DynResidue::new(&g1.0 .0, residue_params.clone()),
                DynResidue::new(&g2.0 .0, residue_params.clone()),
            );

            let (res_n1, res_n2) = (
                DynResidue::new(&n1.0 .0, residue_params.clone()),
                DynResidue::new(&n2.0 .0, residue_params.clone()),
            );

            let res = res_g1 * i1 * res_n2 + res_g2 * i2 * res_n1;

            HenselCode {
                params: residue_params,
                res,
            }
        }
    }
    /// add two HenselCodes
    impl<const L: usize> Add<HenselCode<L>> for HenselCode<L> {
        type Output = HenselCode<L>;
        fn add(self, other: HenselCode<L>) -> HenselCode<L> {
            &self + &other
        }
    }
    /// multiply two HenselCodes
    impl<const L: usize> Mul<HenselCode<L>> for HenselCode<L> {
        type Output = HenselCode<L>;
        fn mul(self, other: HenselCode<L>) -> HenselCode<L> {
            &self * &other
        }
    }

    /// add two &HenselCodes
    impl<'a, 'b, const L: usize> Add<&'b HenselCode<L>> for &'a HenselCode<L> {
        type Output = HenselCode<L>;
        fn add(self, other: &'b HenselCode<L>) -> HenselCode<L> {
            if self.modulus() != other.modulus() {
                panic!("cannot add '{}' and '{}'", self, other);
            }
            HenselCode {
                params: self.params,
                res: self.res + other.res,
            }
        }
    }
    /// multiply two &HenselCodes
    impl<'a, 'b, const L: usize> Mul<&'b HenselCode<L>> for &'a HenselCode<L> {
        type Output = HenselCode<L>;
        fn mul(self, other: &'b HenselCode<L>) -> HenselCode<L> {
            if self.modulus() != other.modulus() {
                panic!("cannot add '{}' and '{}'", self, other);
            }
            HenselCode {
                params: self.params,
                res: self.res * other.res,
            }
        }
    }

    /// pretty-print HenselCode
    impl<const L: usize> fmt::Display for HenselCode<L> {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{} (mod {})", self.to_big_int(), self.modulus())
        }
    }

    /// given a prime `p` and a rational `r = num/denom`, where p does not divide denom,
    /// returns `r (mod p)`
    impl<const L: usize> From<(BigInt<L>, Rational<L>)> for HenselCode<L> {
        fn from(params: (BigInt<L>, Rational<L>)) -> Self {
            let (g, r) = params;
            let params = DynResidueParams::new(&g.0 .0);
            let denom = DynResidue::<L>::new(&r.denom.0 .0, params);
            let num = DynResidue::<L>::new(&r.num.0 .0, params);
            let (id, _) = denom.invert();
            let res = id * num;
            HenselCode { params, res }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::hensel_code::new_hensel_code;

    use super::{big_int::BigInt, hensel_code::HenselCode, rational::Rational, *};

    // #[test]
    // fn finds_modular_inverses() {
    //     let inv = modular_inverses(BigInt::from(37), BigInt::from(5));
    //     assert_eq!(inv, (-2, 15)); // -2*37 + 15*5 = 1
    // }

    // #[test]
    // fn finds_modular_inverses_big() {
    //     let inv = modular_inverses(374533, 52343);
    //     assert_eq!(inv, (-19870, 142177));
    // }

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
            assert_eq!(hc.to_big_int().0 .0, (id_hc * &n_hc).to_big_int().0 .0);
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
            assert_eq!(hc.to_big_int().0 .0, (id_hc * &n_hc).to_big_int().0 .0);
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
