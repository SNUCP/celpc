# armortized-pok

Rust Code for [HSS23].

## How to execute

You need to install nightly toolchain of rustc.
```
$ RUSTFLAGS="-C -target-cpu=native" cargo run --release
keygen time: 3.933953ms
commit time: 234.62µs
open time: 1.975185ms
commit ok? true
ajtai prove time: 1.851098356s
ajtai verify time: 137.031944ms
ajtai verify ok? true
r1cs prove time: 35.89997252s
r1cs verify time: 2.205492363s
r1cs verify ok? true
```
