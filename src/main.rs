use celpc::prelude::*;
use rug::Integer;
use std::time::*;

fn main() {
    let params = Parameters::new(ParametersLiteral::logn_19_logp_256());
    let ck = AjtaiCommitKey::new(&params);
    let mut prover = Prover::new(&params, &ck);
    let mut verifier = Verifier::new(&params, &ck);

    let mut v = vec![Integer::ZERO; params.degree()];
    for i in 0..params.degree() {
        prover.oracle.read_big_mod_assign(&mut v[i]);
    }
    let mut x = Integer::ZERO;
    prover.oracle.read_big_mod_assign(&mut x);

    let now = Instant::now();
    let (com, open) = prover.commit(&v);
    let open_pf = prover.prove_opening(&com, &open);
    let eval_pf = prover.evaluate(x.clone(), &open);
    println!("commit time: {:?}", now.elapsed());

    let now = Instant::now();
    let (com_plain, open_plain) = prover.commit_plain(&v);
    let open_pf_plain = prover.prove_opening_plain(&com_plain, &open_plain);
    let eval_pf_plain = prover.evaluate_plain(x.clone(), &open_plain);
    println!("commit_plain time: {:?}", now.elapsed());

    let now = Instant::now();
    let open_pf_vf = verifier.verify_opening_proof(&com, &open_pf);
    println!("verify_opening_proof time: {:?}", now.elapsed());
    println!("open_pf_vf: {:?}", open_pf_vf);

    let now = Instant::now();
    let eval_pf_vf = verifier.verify_evaluation_proof(x.clone(), &com, &eval_pf);
    println!("verify_evaluation_proof time: {:?}", now.elapsed());
    println!("eval_pf_vf: {:?}", eval_pf_vf);

    let now = Instant::now();
    let open_pf_plain_vf = verifier.verify_opening_proof_plain(&com_plain, &open_pf_plain);
    println!("verify_opening_proof_plain time: {:?}", now.elapsed());
    println!("open_pf_plain_vf: {:?}", open_pf_plain_vf);

    let now = Instant::now();
    let eval_pf_plain_vf =
        verifier.verify_evaluation_proof_plain(x.clone(), &com_plain, &eval_pf_plain);
    println!("verify_evaluation_proof_plain time: {:?}", now.elapsed());
    println!("eval_pf_plain_vf: {:?}", eval_pf_plain_vf);
}
