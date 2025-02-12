#![allow(non_snake_case)]

use std::time::Instant;

use ethnum::U256;
use polycom::{CommitKey, Parameters, PolynomialProver, PolynomialVerifier};

fn main() {
    with_params(Parameters::N_19());
    with_params(Parameters::N_21());
    with_params(Parameters::N_23());
    with_params(Parameters::N_25());
}

fn with_params(pp: Parameters) {
    println!(
        "Current Parameters: N = 2^{} with p = {}^{} + 1.",
        pp.N.ilog2(),
        pp.b,
        pp.kap
    );

    let now = Instant::now();
    let ck = CommitKey::new(&pp, &[0, 0, 0, 0]);
    let mut prover = PolynomialProver::new(&pp, &ck);
    let mut verifier = PolynomialVerifier::new(&pp, &ck);
    println!("setup: {:?}", now.elapsed());

    let p = vec![U256::ONE; pp.N];
    let x = U256::ONE;

    let now = Instant::now();
    let pc = prover.commit(&p);
    println!("commit: {:?}", now.elapsed());

    let now = Instant::now();
    let _ = prover.commit_nozk(&p);
    println!("commit (no zk): {:?}", now.elapsed());

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
    let _ = prover.prove_nozk(&pc);
    println!("open (no zk): {:?}", now.elapsed());

    let now = Instant::now();
    let v = verifier.verify(&pc, &op);
    println!("open_verify: {:?}", now.elapsed());
    println!("open_verify result: {}", v);
}
