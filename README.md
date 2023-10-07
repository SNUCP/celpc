# armortized-pok

Rust Code for [HSS23].

## How to execute

You need to install nightly toolchain of rustc.
```
$ RUSTFLAGS="-C -target-cpu=native" cargo run --release
keygen time: 3.998551ms
commit time: 217.079µs
open time: 1.838292ms
commit ok? true
ajtai prove time: 1.447710452s
ajtai verify time: 114.600949ms
ajtai verify ok? true
r1cs prove time: 35.526223249s
r1cs verify time: 1.033504665s
r1cs verify ok? true
```
