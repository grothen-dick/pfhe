use crypto_bigint::{Checked, Wrapping, U256};

/// default size of a BigInt (in LIMBS)
pub const DEFAULT_LIMBS: usize = U256::LIMBS;

/// an enum to specify wether to use Wrapping of Checked arithmetics
pub enum Arithmetics {
    Wrapping,
    Checked,
}

/// the type of arithmetics used globally
pub const ARITHMETICS: Arithmetics = Arithmetics::Wrapping;

/// a macro to create a new BigInt with arithmetics specified by ARITHMETICS
#[macro_export]
macro_rules! new_usable_bigint {
    ($expression: expr, $arithmetics: path) => {
        $arithmetics($expression)
    };
}

/// A trait to retrieve the current size of a const-generic struct
pub trait Bounded {
    const L: usize;

    fn size(&self) -> usize {
        Self::L
    }
}
