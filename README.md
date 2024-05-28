# CELPC: Concretely Efficient Lattice-based Polynomial Commitment

Rust Code for [[HSS24](https://eprint.iacr.org/2024/306)] "Concretely Efficient Lattice-based Polynomial Commitment".

You need to install nightly toolchain of rustc to run the code.
```
$ RUSTFLAGS="-C target-cpu=native" cargo run --release
Current Parameters: N = 2^19 with p = 63388^16 + 1.
commit: 3.672398851s
commit (no zk): 502.869079ms
evaluate: 37.417646ms
evaluate result: 524288 ?= 524288
eval_verify: 17.57424ms
eval_verify result: true
open: 820.32421ms
open (no zk): 410.326699ms
open_verify: 121.148132ms
open_verify result: true
...
```
