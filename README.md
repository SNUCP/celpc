# CELPC: Concretely Efficient Lattice-based Polynomial Commitment

Rust Code for [[HSS24](https://eprint.iacr.org/2024/306)] "Concretely Efficient Lattice-based Polynomial Commitment".

You need to install nightly toolchain of rustc to run the code.
```
$ RUSTFLAGS="-C target-cpu=native" cargo run --release
Current Parameters: N = 2^19 with p = 63388^16 + 1.
setup: 1.433771978s
commit: 2.656603899s
commit (no zk): 501.154737ms
evaluate: 38.001503ms
evaluate result: 524288 ?= 524288
eval_verify: 17.819188ms
eval_verify result: true
open: 802.689522ms
open (no zk): 406.850866ms
open_verify: 125.980332ms
open_verify result: true
...
```
