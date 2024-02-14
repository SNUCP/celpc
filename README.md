# ELSA-PC: Concretely Efficient Lattice-based Polynomial Commitment

Rust Code for [HSS24] "Concretely Efficient Lattice-based Polynomial Commitment".

You need to install nightly toolchain of rustc to run the code. Currently works only in Rust Nightly 1.77.0.
```
$ RUSTFLAGS="-C target-cpu=native" cargo run --release
Current Parameters: N = 2^19 with p = 63388^16 + 1.
commit: 4.187749949s
commit (no zk): 507.360968ms
evaluate: 37.907132ms
evaluate result: 524288 ?= 524288
eval_verify: 17.685859ms
eval_verify result: true
open: 738.555925ms
open (no zk): 420.328664ms
open_verify: 122.901533ms
open_verify result: true
...
```
