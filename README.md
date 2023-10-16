# ELSA: Efficient Lattice-based Sublinear Arguments for R1CS without Aborts

Rust Code for [HSS23] "Efficient Lattice-based Sublinear Arguments for R1CS without Aborts".

You need to install nightly toolchain of rustc to run the code.
```
$ RUSTFLAGS="-C target-cpu=native" cargo run --release
keygen time: 3.915143ms
commit time: 218.024µs
open time: 1.963427ms
commit ok? true
ajtai prove time: 1.498060698s
ajtai verify time: 114.837694ms
ajtai verify ok? true
r1cs prove time: 33.257246645s
r1cs verify time: 1.030204535s
r1cs verify ok? true
```
