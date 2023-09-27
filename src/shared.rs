use crypto_bigint::U256;

/// Default size of a BigInt (in LIMBS)
pub const DEFAULT_LIMBS: usize = U256::LIMBS;

/// Trait to retrieve the current size of a const-generic struct
pub trait Bounded {
    const L: usize;

    fn size(&self) -> usize {
        Self::L
    }
}
