#[macro_export]
macro_rules! impl_arithmetic_op {
    ($trait: ident, $function: ident, $op: tt) => {
        impl<'a, 'b, const L: usize> $trait<&'b BigInt<L>> for &'a BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: &'b BigInt<L>) -> BigInt<L> {
                BigInt(self.0 $op other.0)
            }
        }

        impl<const L: usize> $trait<BigInt<L>> for BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: BigInt<L>) -> BigInt<L> {
                &self $op &other
            }
        }

    };
}

#[macro_export]
macro_rules! impl_divlike_op {
    ($trait: ident, $function: ident, $op: tt) => {
        impl<'a, 'b, const L: usize> $trait<&'b BigInt<L>> for &'a BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: &'b BigInt<L>) -> BigInt<L> {
                BigInt(self.0 $op NonZero::new(other.to_uint()).unwrap())
            }
        }

        impl<const L: usize> $trait<BigInt<L>> for BigInt<L> {
            type Output = BigInt<L>;
            fn $function(self, other: BigInt<L>) -> BigInt<L> {
                &self $op &other
            }
        }

    };
}
