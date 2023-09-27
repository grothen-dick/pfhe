# :closed_lock_with_key: P-adic Fully Homomorphic Encryption

![main branch ci status](https://github.com/p-adic-fhe/pfhe/actions/workflows/ci.yml/badge.svg?branch=main)

## ⚠️ WARNING! ⚠️

This repo contains highly experimental code.
Expect rapid iteration.

## :eyeglasses: Description

This repository is a library which provides: 

- methods for a FHE (Fully Homomorphic Encryption) scheme that uses p-adic numbers via Hensel codes
- conversion of rational numbers to and from Hensel codes 

Our method is based on this [research paper](https://eprint.iacr.org/2021/1281.pdf) by David W. H. A. da Silva, Luke Harmon, Gaetan Delavignette, and Carlos Araujo.

## :hammer: Build

```rust
cargo build
```

## :100: Test

```rust
cargo test
```
