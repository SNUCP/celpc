#![allow(non_snake_case)]

use std::time::Instant;

use apok::csprng::*;
use apok::*;
use ethnum::U256;

fn main() {
    let params = Parameters::default();
    let mut ecd = Encoder::new(&params);
    let mut us = UniformSampler::new();
    let mut gs = KarneySampler::new();

    let mut m = vec![params.ringq.new_poly(); params.l];
    let mut msg = vec![U256::ZERO; params.m];
    for i in 0..params.l {
        for j in 0..params.m {
            msg[j] = us.sample_u256() % params.p;
        }
        ecd.encode_randomized_assign(&msg, params.s1, &mut m[i]);
    }

    let mut mu = vec![params.ringq.new_poly(); params.munu];
    for i in 0..params.munu {
        gs.sample_poly_assign(&params.ringq, 0.0, params.sig1, &mut mu[i]);
    }

    let now = Instant::now();
    let ck = CommitKey::new(&params, &[0u8; 32]);
    println!("keygen time: {:?}", now.elapsed());

    let commiter = Comitter::new(&params, &ck);

    let now = Instant::now();
    let commit = commiter.commit(&m, &mu);
    println!("commit time: {:?}", now.elapsed());

    let now = Instant::now();
    let ok = commiter.open(&m, &mu, (f64::exp2(20.0), f64::exp(20.0)), &commit);
    println!("open time: {:?}", now.elapsed());
    println!("commit ok? {:?}", ok);

    let mut ajprover = AjtaiProver::new(&params, &ck);

    let ajm = vec![m.clone(); 4 * params.k + 2];
    let ajmu = vec![mu.clone(); 4 * params.k + 2];
    let ajc = vec![commit.clone(); 4 * params.k + 2];

    let now = Instant::now();
    let proof = ajprover.prove(&ajm, &ajmu, &ajc);
    println!("ajtai prove time: {:?}", now.elapsed());

    let mut ajverifier = AjtaiVerifier::new(&params, &ck);

    let now = Instant::now();
    let ok = ajverifier.verify(&ajc, &proof);
    println!("ajtai verify time: {:?}", now.elapsed());
    println!("ajtai verify ok? {:?}", ok);

    let A = SparseMatrix::new_identity(params.M, params.p);
    let B = SparseMatrix::new_identity(params.M, params.p);
    let C = SparseMatrix::new_identity(params.M, params.p);

    let a = vec![U256::from(1u32); params.M];
    let b = vec![U256::from(2u32); params.M];
    let c = vec![U256::from(5u32); params.M];

    let tu = vec![U256::from(1u32); params.m];
    let mut t = vec![vec![params.ringq.new_ntt_poly(); params.l]; params.k];
    let mut tau = vec![vec![params.ringq.new_ntt_poly(); params.munu]; params.k];
    let mut tcommit = vec![vec![params.ringq.new_ntt_poly(); params.mu]; params.k];
    for i in 0..params.k {
        for j in 0..params.l {
            ecd.encode_randomized_assign(&tu, params.s1, &mut t[i][j]);
        }
        for j in 0..params.munu {
            gs.sample_poly_assign(&params.ringq, 0.0, params.sig1, &mut tau[i][j]);
        }
        commiter.commit_assign(&t[i], &tau[i], &mut tcommit[i]);
    }

    let mut r1csprover = R1CSProver::new(&params, &ck);

    let now = Instant::now();
    let proof = r1csprover.prove(&A, &B, &C, &a, &b, &c, &t, &tau, &tcommit);
    println!("r1cs prove time: {:?}", now.elapsed());

    let mut r1csverifier = R1CSVerifier::new(&params, &ck);

    let now = Instant::now();
    let ok = r1csverifier.verify(&A, &B, &C, &a, &b, &c, &tcommit, &proof);
    println!("r1cs verify time: {:?}", now.elapsed());
    println!("r1cs verify ok? {:?}", ok);
}
