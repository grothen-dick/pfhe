use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pfhe::{
    bigint::{BigIntTrait, WrappingCryptoBigInt},
    crypto_parameters::CryptographicParameters,
    hensel_code::HenselCode,
    rational::Rational,
};
use std::time::Instant;

pub fn criterion_benchmark(c: &mut Criterion) {
    type T = WrappingCryptoBigInt;
    let p = T::from_u128(7919); // thanks wikipedia for this prime

    let r = Rational::<T> {
        num: T::from_u128(6),
        denom: T::from_u128(1),
    };

    c.bench_function("encode rational to hensel code", |b| {
        b.iter(|| HenselCode::from(black_box((&p, &r))))
    });

    const L1: usize = 170;
    type T1 = WrappingCryptoBigInt<L1>;

    println!("generating crypto params ...");
    let now = Instant::now();
    let crypto_params = CryptographicParameters::<T1>::from_params(8, 3);
    println!("crypto params generated in {:.2?}", now.elapsed());

    let message: Rational<T1> = Rational {
        num: T1::from_u128(43),
        denom: T1::from_u128(44),
    };

    c.bench_function("encrypt rational using large prime number", |b| {
        b.iter(|| crypto_params.encrypt(black_box(message.clone())))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
