use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pfhe::{
    bigint::BigInt, crypto_parameters::CryptographicParameters, hensel_code::HenselCode,
    rational::Rational, shared::DEFAULT_LIMBS,
};
use std::time::Instant;

pub fn criterion_benchmark(c: &mut Criterion) {
    const L: usize = DEFAULT_LIMBS;
    let p: BigInt<L> = BigInt::<L>::from(7919 as u128); // thanks wikipedia for this prime

    let r = Rational::<L> {
        num: BigInt::<L>::from(6 as u128),
        denom: BigInt::<L>::from(1 as u128),
    };

    c.bench_function("encode rational to hensel code", |b| {
        b.iter(|| HenselCode::from(black_box((&p, &r))))
    });

    const L1: usize = 40;

    println!("generating crypto params ...");
    let now = Instant::now();
    let crypto_params = CryptographicParameters::<L1>::from_params(4, 3);
    println!("crypto params generated in {:.2?}", now.elapsed());

    let message: Rational<L1> = Rational {
        num: BigInt::<L1>::from(43),
        denom: BigInt::<L1>::from(44),
    };

    c.bench_function("encrypt rational using large prime number", |b| {
        b.iter(|| crypto_params.encrypt(black_box(message.clone())))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
