use crypto_bigint::U256;
/// default size of a BigInt (in LIMBS)
pub const DEFAULT_LIMBS: usize = U256::LIMBS;

/// A trait to retrieve the current size of a const-generic struct
pub trait Bounded {
    const L: usize;

    fn size(&self) -> usize {
        Self::L
    }
}
