use celpc::prelude::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rug::Integer;

fn celpc_benchmark(c: &mut Criterion) {
    let params_literal_list = [
        ParametersLiteral::logn_19_logp_256(),
        ParametersLiteral::logn_21_logp_256(),
        ParametersLiteral::logn_23_logp_256(),
        ParametersLiteral::logn_25_logp_256(),
    ];

    for params_literal in params_literal_list.iter() {
        let params = Parameters::new(params_literal.clone());
        let ck = AjtaiCommitKey::new(&params);
        let mut prover = Prover::new(&params, &ck);
        let mut verifier = Verifier::new(&params, &ck);

        let mut v = vec![Integer::ZERO; params.degree()];
        for i in 0..params.degree() {
            prover.oracle.read_big_mod_assign(&mut v[i]);
        }
        let mut x = Integer::ZERO;
        prover.oracle.read_big_mod_assign(&mut x);

        let mut com = Commitment::new(&params, params.degree());
        let mut open = Opening::new(&params, params.degree());
        let mut open_pf = OpeningProof::new(&params);
        let mut eval_pf = EvaluationProof::new(&params);

        c.bench_function(
            &format!("ParamsLog{}P256_commit", params.degree().ilog2()),
            |b| {
                b.iter(|| {
                    prover.commit_assign(black_box(&v), black_box(&mut com), black_box(&mut open))
                })
            },
        );

        c.bench_function(
            &format!("ParamsLog{}P256_prove_opening", params.degree().ilog2()),
            |b| {
                b.iter(|| {
                    prover.prove_opening_assign(
                        black_box(&com),
                        black_box(&open),
                        black_box(&mut open_pf),
                    )
                })
            },
        );

        c.bench_function(
            &format!("ParamsLog{}P256_evaluate", params.degree().ilog2()),
            |b| {
                b.iter(|| {
                    prover.evaluate_assign(
                        black_box(x.clone()),
                        black_box(&open),
                        black_box(&mut eval_pf),
                    )
                })
            },
        );

        c.bench_function(
            &format!("ParamsLog{}P256_commit_plain", params.degree().ilog2()),
            |b| {
                b.iter(|| {
                    prover.commit_plain_assign(
                        black_box(&v),
                        black_box(&mut com),
                        black_box(&mut open),
                    )
                })
            },
        );

        c.bench_function(
            &format!(
                "ParamsLog{}P256_prove_opening_plain",
                params.degree().ilog2()
            ),
            |b| {
                b.iter(|| {
                    prover.prove_opening_plain_assign(
                        black_box(&com),
                        black_box(&open),
                        black_box(&mut open_pf),
                    )
                })
            },
        );

        c.bench_function(
            &format!("ParamsLog{}P256_evaluate_plain", params.degree().ilog2()),
            |b| {
                b.iter(|| {
                    prover.evaluate_plain_assign(
                        black_box(x.clone()),
                        black_box(&open),
                        black_box(&mut eval_pf),
                    )
                })
            },
        );

        c.bench_function(
            &format!(
                "ParamsLog{}P256_verify_opening_proof",
                params.degree().ilog2()
            ),
            |b| b.iter(|| verifier.verify_opening_proof(black_box(&com), black_box(&open_pf))),
        );

        c.bench_function(
            &format!(
                "ParamsLog{}P256_verify_evaluation_proof",
                params.degree().ilog2()
            ),
            |b| {
                b.iter(|| {
                    verifier.verify_evaluation_proof(
                        black_box(x.clone()),
                        black_box(&com),
                        black_box(&eval_pf),
                    )
                })
            },
        );
    }
}

criterion_group!(benches, celpc_benchmark);
criterion_main!(benches);
