use crate::{ajtai::*, csprng::*, encoder::*, entities::*, parameters::Parameters};
use rug::{Assign, Integer};

/// Verifier is a verifier for polynomial commitments.
pub struct Verifier {
    pub params: Parameters,

    pub oracle: UniformSampler,

    pub encoder: Encoder,
    pub commiter: AjtaiCommiter,
}

impl Verifier {
    /// Creates a new verifier.
    pub fn new(params: &Parameters, ck: &AjtaiCommitKey) -> Verifier {
        Verifier {
            params: params.clone(),

            oracle: UniformSampler::new_oracle(params),
            encoder: Encoder::new(params),
            commiter: AjtaiCommiter::new(params, ck),
        }
    }

    /// Verifies the proof of opening knowledge.
    pub fn verify_opening_proof(&mut self, com: &Commitment, open_pf: &OpeningProof) -> bool {
        let batch_count = com.value.len() - 1;

        self.oracle.reset();
        for i in 0..batch_count {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&com.value[i].value[j]);
            }
        }

        for i in 0..self.params.repetition() {
            let mut norm = 0.0;
            for j in 0..self.params.poly_commit_size() {
                norm += self.params.ringq().norm_sq(&open_pf.response_mask[i][j]);
            }

            for j in 0..self.params.ajtai_rand_size() {
                norm += self.params.ringq().norm_sq(&open_pf.response_rand[i][j]);
            }

            if norm.sqrt() > self.params.open_proof_bound() {
                return false;
            }
        }

        for i in 0..self.params.repetition() {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&open_pf.commitment[i].value[j]);
            }
        }

        let mut challenge =
            vec![vec![self.params.ringq().new_poly(); batch_count]; self.params.repetition()];
        for i in 0..self.params.repetition() {
            for j in 0..batch_count {
                self.oracle.read_monomial_assign(&mut challenge[i][j]);
            }
        }

        let mut commit_response = AjtaiCommitment::new(&self.params);
        let mut commit_combine = AjtaiCommitment::new(&self.params);
        for i in 0..self.params.repetition() {
            self.commiter.commit_assign(
                &open_pf.response_mask[i],
                &open_pf.response_rand[i],
                &mut commit_response,
            );

            commit_combine.clone_from(&open_pf.commitment[i]);
            for j in 0..batch_count {
                for k in 0..self.params.ajtai_size() {
                    self.params.ringq().mul_add_assign(
                        &challenge[i][j],
                        &com.value[j].value[k],
                        &mut commit_combine.value[k],
                    );
                }
            }

            if !commit_response.equals(&commit_combine) {
                return false;
            }
        }

        return true;
    }

    /// Verifies the proof of opening knowledge without hiding.
    pub fn verify_opening_proof_plain(&mut self, com: &Commitment, open_pf: &OpeningProof) -> bool {
        let batch_count = com.value.len() - 2;

        self.oracle.reset();
        for i in 0..batch_count {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&com.value[i].value[j]);
            }
        }

        for i in 0..self.params.repetition() {
            let mut norm = 0.0;
            for j in 0..self.params.poly_commit_size() {
                norm += self.params.ringq().norm_sq(&open_pf.response_mask[i][j]);
            }

            if norm.sqrt() > self.params.open_proof_bound() {
                return false;
            }
        }

        for i in 0..self.params.repetition() {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&open_pf.commitment[i].value[j]);
            }
        }

        let mut challenge =
            vec![vec![self.params.ringq().new_poly(); batch_count]; self.params.repetition()];
        for i in 0..self.params.repetition() {
            for j in 0..batch_count {
                self.oracle.read_monomial_assign(&mut challenge[i][j]);
            }
        }

        let mut commit_response = AjtaiCommitment::new(&self.params);
        let mut commit_combine = AjtaiCommitment::new(&self.params);
        for i in 0..self.params.repetition() {
            self.commiter
                .commit_plain_assign(&open_pf.response_mask[i], &mut commit_response);

            commit_combine.clone_from(&open_pf.commitment[i]);
            for j in 0..batch_count {
                for k in 0..self.params.ajtai_size() {
                    self.params.ringq().mul_add_assign(
                        &challenge[i][j],
                        &com.value[j].value[k],
                        &mut commit_combine.value[k],
                    );
                }
            }

            if !commit_response.equals(&commit_combine) {
                return false;
            }
        }

        return true;
    }

    /// Verifies the proof of evaluation.
    pub fn verify_evaluation_proof(
        &mut self,
        x: Integer,
        com: &Commitment,
        eval_pf: &EvaluationProof,
    ) -> bool {
        let mut norm = 0.0;
        for i in 0..self.params.poly_commit_size() {
            norm += self.params.ringq().norm_sq(&eval_pf.mask[i]);
        }

        for i in 0..self.params.ajtai_rand_size() {
            norm += self.params.ringq().norm_sq(&eval_pf.rand[i]);
        }

        if norm.sqrt() > self.params.eval_bound() {
            return false;
        }

        let commit_count = com.value.len() - 2;

        let x_ecd = self
            .encoder
            .encode(&[x.clone().modulo(&self.params.modulus())]);
        let x_pow_n = x
            .clone()
            .pow_mod(
                &Integer::from(self.params.bigint_commit_size()),
                self.params.modulus(),
            )
            .unwrap();
        let mut x_pow_buf = Integer::from(1);
        let mut x_pow_ecd = vec![self.params.ringq().new_poly(); commit_count];
        for i in 0..commit_count {
            self.encoder
                .encode_assign(&[x_pow_buf.clone()], &mut x_pow_ecd[i]);
            x_pow_buf *= &x_pow_n;
            x_pow_buf %= self.params.modulus();
        }

        let commit_eval_pf = self.commiter.commit(&eval_pf.mask, &eval_pf.rand);
        let mut commit_combine = AjtaiCommitment::new(&self.params);
        commit_combine.clone_from(&com.value[commit_count + 1]);
        for j in 0..self.params.ajtai_size() {
            self.params.ringq().mul_add_assign(
                &x_ecd,
                &com.value[commit_count].value[j],
                &mut commit_combine.value[j],
            );
        }
        for i in 0..commit_count {
            for j in 0..self.params.ajtai_size() {
                self.params.ringq().mul_add_assign(
                    &x_pow_ecd[i],
                    &com.value[i].value[j],
                    &mut commit_combine.value[j],
                );
            }
        }

        if !commit_eval_pf.equals(&commit_combine) {
            return false;
        }

        let mut x_pow_basis = vec![Integer::ZERO; self.params.bigint_commit_size()];
        x_pow_buf.assign(1);
        for i in 0..self.params.bigint_commit_size() {
            x_pow_basis[i].assign(&x_pow_buf);
            x_pow_buf *= &x;
            x_pow_buf %= self.params.modulus();
        }

        let mask_dcd = self.encoder.decode_chunk(&eval_pf.mask);
        let mut value_combine = Integer::ZERO;
        let mut prod_buf = Integer::ZERO;
        for i in 0..self.params.bigint_commit_size() {
            prod_buf.assign(&mask_dcd[i]);
            prod_buf *= &x_pow_basis[i];
            value_combine += &prod_buf;
        }
        value_combine %= self.params.modulus();

        return eval_pf.value == value_combine;
    }

    /// Verifies the proof of evaluation without hiding.
    pub fn verify_evaluation_proof_plain(
        &mut self,
        x: Integer,
        com: &Commitment,
        eval_pf: &EvaluationProof,
    ) -> bool {
        let mut norm = 0.0;
        for i in 0..self.params.poly_commit_size() {
            norm += self.params.ringq().norm_sq(&eval_pf.mask[i]);
        }

        if norm.sqrt() > self.params.eval_bound() {
            return false;
        }

        let commit_count = com.value.len() - 2;

        let x_pow_n = x
            .clone()
            .pow_mod(
                &Integer::from(self.params.bigint_commit_size()),
                self.params.modulus(),
            )
            .unwrap();
        let mut x_pow_buf = Integer::from(1);
        let mut x_pow_ecd = vec![self.params.ringq().new_poly(); commit_count];
        for i in 0..commit_count {
            self.encoder
                .encode_assign(&[x_pow_buf.clone()], &mut x_pow_ecd[i]);
            x_pow_buf *= &x_pow_n;
            x_pow_buf %= self.params.modulus();
        }

        let commit_eval_pf = self.commiter.commit(&eval_pf.mask, &eval_pf.rand);
        let mut commit_combine = AjtaiCommitment::new(&self.params);
        for j in 0..self.params.ajtai_size() {
            self.params.ringq().mul_assign(
                &x_pow_ecd[0],
                &com.value[0].value[j],
                &mut commit_combine.value[j],
            );
        }
        for i in 1..commit_count {
            for j in 0..self.params.ajtai_size() {
                self.params.ringq().mul_add_assign(
                    &x_pow_ecd[i],
                    &com.value[i].value[j],
                    &mut commit_combine.value[j],
                );
            }
        }

        if !commit_eval_pf.equals(&commit_combine) {
            return false;
        }

        let mut x_pow_basis = vec![Integer::ZERO; self.params.bigint_commit_size()];
        x_pow_buf.assign(1);
        for i in 0..self.params.bigint_commit_size() {
            x_pow_basis[i].assign(&x_pow_buf);
            x_pow_buf *= &x;
            x_pow_buf %= self.params.modulus();
        }

        let mask_dcd = self.encoder.decode_chunk(&eval_pf.mask);
        let mut value_combine = Integer::ZERO;
        let mut prod_buf = Integer::ZERO;
        for i in 0..self.params.bigint_commit_size() {
            prod_buf.assign(&mask_dcd[i]);
            prod_buf *= &x_pow_basis[i];
            value_combine += &prod_buf;
        }
        value_combine %= self.params.modulus();

        return eval_pf.value == value_combine;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use rug::Integer;

    #[test]
    fn test_celpc() {
        let params = Parameters::new(ParametersLiteral::logn_19_logp_256());
        let ck = AjtaiCommitKey::new(&params);
        let mut prover = Prover::new(&params, &ck);
        let mut verifier = Verifier::new(&params, &ck);
        let mut oracle = UniformSampler::new_oracle(&params);

        let mut v = vec![Integer::ZERO; params.bigint_commit_size()];
        for i in 0..v.len() {
            oracle.read_big_mod_assign(&mut v[i]);
        }
        let mut x = Integer::ZERO;
        oracle.read_big_mod_assign(&mut x);

        let (com, open) = prover.commit(&v);
        let open_pf = prover.prove_opening(&com, &open);
        let eval_pf = prover.evaluate(x.clone(), &open);

        assert!(verifier.verify_opening_proof(&com, &open_pf));
        assert!(verifier.verify_evaluation_proof(x.clone(), &com, &eval_pf));

        let mut value_real = Integer::ZERO;
        for i in (0..v.len()).rev() {
            value_real *= &x;
            value_real += &v[i];
            value_real %= params.modulus();
        }

        assert_eq!(eval_pf.value, value_real);
    }

    #[test]
    fn test_celpc_plain() {
        let params = Parameters::new(ParametersLiteral::logn_19_logp_256());
        let ck = AjtaiCommitKey::new(&params);
        let mut prover = Prover::new(&params, &ck);
        let mut verifier = Verifier::new(&params, &ck);
        let mut oracle = UniformSampler::new_oracle(&params);

        let mut v = vec![Integer::ZERO; params.bigint_commit_size()];
        for i in 0..v.len() {
            oracle.read_big_mod_assign(&mut v[i]);
        }
        let mut x = Integer::ZERO;
        oracle.read_big_mod_assign(&mut x);

        let (com, open) = prover.commit_plain(&v);
        let open_pf = prover.prove_opening_plain(&com, &open);
        let eval_pf = prover.evaluate_plain(x.clone(), &open);

        assert!(verifier.verify_opening_proof_plain(&com, &open_pf));
        assert!(verifier.verify_evaluation_proof_plain(x.clone(), &com, &eval_pf));

        let mut value_real = Integer::ZERO;
        for i in (0..v.len()).rev() {
            value_real *= &x;
            value_real += &v[i];
            value_real %= params.modulus();
        }

        assert_eq!(eval_pf.value, value_real);
    }
}
