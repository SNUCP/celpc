#![allow(non_snake_case)]

use std::time::Instant;

use apok::{CommitKey, Parameters, PolynomialProver, PolynomialVerifier};
use ethnum::U256;

fn main() {
    let pp = Parameters::N_19();

    let ck = CommitKey::new(&pp, &[0, 0, 0, 0]);
    let mut prover = PolynomialProver::new(&pp, &ck);
    let mut verifier = PolynomialVerifier::new(&pp, &ck);

    let p = vec![U256::ONE; pp.N];
    let x = U256::ONE;

    let now = Instant::now();
    let pc = prover.commit(&p);
    println!("commit: {:?}", now.elapsed());

    let now = Instant::now();
    let (y, ep) = prover.evaluate(x, &pc);
    println!("evaluate: {:?}", now.elapsed());
    println!("evaluate result: {:?} ?= {:?}", y, pp.N);

    let now = Instant::now();
    let v = verifier.verify_evaluation(x, y, &pc, &ep);
    println!("eval_verify: {:?}", now.elapsed());
    println!("eval_verify result: {}", v);

    let now = Instant::now();
    let op = prover.prove(&pc);
    println!("open: {:?}", now.elapsed());

    let now = Instant::now();
    let v = verifier.verify(&pc, &op);
    println!("open_verify: {:?}", now.elapsed());
    println!("open_verify result: {}", v);
}
