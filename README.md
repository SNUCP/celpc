# CELPC: Concretely Efficient Lattice-based Polynomial Commitment

Rust Code for [[HSS24](https://eprint.iacr.org/2024/306)] "Concretely Efficient Lattice-based Polynomial Commitment".

You need to install nightly toolchain of rustc to run the code.
```
$ RUSTFLAGS="-C target-cpu=native" cargo run --release
Current Parameters: N = 2^19 with p = 63388^16 + 1.
setup: 1.725607199s
commit: 1.287635716s
commit (no zk): 516.111519ms
evaluate: 38.460491ms
evaluate result: 524288 ?= 524288
eval_verify: 18.133062ms
eval_verify result: true
open: 778.639485ms
open (no zk): 412.582167ms
open_verify: 123.850852ms
open_verify result: true
...
```
