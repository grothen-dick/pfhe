use criterion::{black_box, criterion_group, criterion_main, Criterion};
use pfhe::{bigint::BigInt, hensel_code::HenselCode, rational::Rational, shared::DEFAULT_LIMBS};

pub fn criterion_benchmark(c: &mut Criterion) {
    const L: usize = DEFAULT_LIMBS;
    let p: BigInt = BigInt::<L>::from(7919 as u128); // thanks wikipedia for this prime

    let r = Rational::<L> {
        num: BigInt::<L>::from(6 as u128),
        denom: BigInt::<L>::from(1 as u128),
    };

    c.bench_function("encode rational", |b| {
        b.iter(|| HenselCode::from(black_box((&p, &r))))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
