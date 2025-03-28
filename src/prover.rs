use crate::{ajtai::*, csprng::*, encoder::*, entities::*, parameters::Parameters};
use rug::{Assign, Integer};

/// Prover is a prover for polynomial commitments.
pub struct Prover {
    pub params: Parameters,

    pub uniform_sampler: UniformSampler,
    pub oracle: UniformSampler,
    pub gaussian_sampler: COSACSampler,
    pub commit_rand_gaussian_sampler: TwinCDTSampler,

    pub encoder: Encoder,
    pub commit_encoder: EncoderFixed,
    pub commiter: AjtaiCommiter,
}

impl Prover {
    /// Creates a new prover.
    pub fn new(params: &Parameters, ck: &AjtaiCommitKey) -> Prover {
        Prover {
            params: params.clone(),

            uniform_sampler: UniformSampler::new(params),
            oracle: UniformSampler::new(params),
            gaussian_sampler: COSACSampler::new(params),
            commit_rand_gaussian_sampler: TwinCDTSampler::new(params, params.commit_rand_std_dev()),

            encoder: Encoder::new(params),
            commit_encoder: EncoderFixed::new(params, params.commit_std_dev()),
            commiter: AjtaiCommiter::new(params, ck),
        }
    }

    /// Commits a polynomial.
    pub fn commit(&mut self, coeffs: &[Integer]) -> (Commitment, Opening) {
        let mut com_out = Commitment::new(&self.params, coeffs.len());
        let mut open_out = Opening::new(&self.params, coeffs.len());
        self.commit_assign(coeffs, &mut com_out, &mut open_out);
        return (com_out, open_out);
    }

    /// Commits a polynomial and writes to com_out and open_out.
    pub fn commit_assign(
        &mut self,
        coeffs: &[Integer],
        com_out: &mut Commitment,
        open_out: &mut Opening,
    ) {
        let commit_count = coeffs.len() / self.params.bigint_commit_size();

        for (i, chunk) in coeffs
            .chunks_exact(self.params.bigint_commit_size())
            .enumerate()
        {
            self.commit_encoder
                .encode_randomize_chunk_assign(chunk, &mut open_out.mask[i]);
            for j in 0..self.params.ajtai_rand_size() {
                self.commit_rand_gaussian_sampler
                    .read_poly_assign(0.0, &mut open_out.rand[i][j]);
            }
            self.commiter.commit_assign(
                &open_out.mask[i],
                &open_out.rand[i],
                &mut com_out.value[i],
            );
        }

        let mut blind0 = vec![Integer::ZERO; self.params.bigint_commit_size()];
        let mut blind1 = vec![Integer::ZERO; self.params.bigint_commit_size()];
        blind0[self.params.bigint_commit_size() - 1].assign(0);
        blind1[0].assign(0);
        for i in 0..self.params.bigint_commit_size() - 1 {
            self.uniform_sampler.read_big_mod_assign(&mut blind0[i]);
            blind1[i + 1] -= &blind0[i];
            blind1[i + 1] += self.params.modulus();
        }

        self.commit_encoder
            .encode_randomize_chunk_assign(&blind0, &mut open_out.mask[commit_count]);
        for j in 0..self.params.ajtai_rand_size() {
            self.commit_rand_gaussian_sampler
                .read_poly_assign(0.0, &mut open_out.rand[commit_count][j]);
        }
        self.commiter.commit_assign(
            &open_out.mask[commit_count],
            &open_out.rand[commit_count],
            &mut com_out.value[commit_count],
        );

        let blind_std_dev = ((commit_count + 2) as f64).sqrt() * self.params.blind_std_dev();
        let blind_rand_std_dev =
            ((commit_count + 2) as f64).sqrt() * self.params.blind_rand_std_dev();
        self.encoder.encode_randomize_chunk_assign(
            &blind1,
            blind_std_dev,
            &mut open_out.mask[commit_count + 1],
        );
        for j in 0..self.params.ajtai_rand_size() {
            self.gaussian_sampler.read_poly_assign(
                0.0,
                blind_rand_std_dev,
                &mut open_out.rand[commit_count + 1][j],
            );
        }
        self.commiter.commit_assign(
            &open_out.mask[commit_count + 1],
            &open_out.rand[commit_count + 1],
            &mut com_out.value[commit_count + 1],
        );
    }

    /// Commits a polynomial without hiding.
    pub fn commit_plain(&mut self, coeffs: &[Integer]) -> (Commitment, Opening) {
        let mut com_out = Commitment::new(&self.params, coeffs.len());
        let mut open_out = Opening::new(&self.params, coeffs.len());
        self.commit_plain_assign(coeffs, &mut com_out, &mut open_out);
        return (com_out, open_out);
    }

    /// Commits a polynomial without hiding and writes to com_out and open_out.
    pub fn commit_plain_assign(
        &mut self,
        coeffs: &[Integer],
        com_out: &mut Commitment,
        open_out: &mut Opening,
    ) {
        for (i, chunk) in coeffs
            .chunks_exact(self.params.bigint_commit_size())
            .enumerate()
        {
            self.encoder
                .encode_chunk_assign(chunk, &mut open_out.mask[i]);
            self.commiter
                .commit_plain_assign(&open_out.mask[i], &mut com_out.value[i]);
        }
    }

    /// Proves the opening.
    pub fn prove_opening(&mut self, com: &Commitment, open: &Opening) -> OpeningProof {
        let mut open_pf_out = OpeningProof::new(&self.params);
        self.prove_opening_assign(com, open, &mut open_pf_out);
        return open_pf_out;
    }

    /// Proves the opening and writes to open_pf_out.
    pub fn prove_opening_assign(
        &mut self,
        com: &Commitment,
        open: &Opening,
        open_pf_out: &mut OpeningProof,
    ) {
        let batch_count = open.mask.len() - 1;

        self.oracle.reset();
        for i in 0..batch_count {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&com.value[i].value[j]);
            }
        }

        let opening_proof_std_dev =
            ((batch_count + 1) as f64).sqrt() * self.params.opening_proof_std_dev();
        let opening_proof_rand_std_dev =
            ((batch_count + 1) as f64).sqrt() * self.params.opening_proof_rand_std_dev();

        let mut rand_mask = vec![Integer::ZERO; self.params.bigint_commit_size()];
        for i in 0..self.params.repetition() {
            for j in 0..self.params.poly_commit_size() {
                self.uniform_sampler.read_big_mod_assign(&mut rand_mask[j]);
            }
            self.encoder.encode_randomize_chunk_assign(
                &rand_mask,
                opening_proof_std_dev,
                &mut open_pf_out.response_mask[i],
            );
            for j in 0..self.params.ajtai_rand_size() {
                self.gaussian_sampler.read_poly_assign(
                    0.0,
                    opening_proof_rand_std_dev,
                    &mut open_pf_out.response_rand[i][j],
                );
            }

            self.commiter.commit_assign(
                &open_pf_out.response_mask[i],
                &open_pf_out.response_rand[i],
                &mut open_pf_out.commitment[i],
            );
        }

        for i in 0..self.params.repetition() {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&open_pf_out.commitment[i].value[j]);
            }
        }

        let mut challenge =
            vec![vec![self.params.ringq().new_poly(); batch_count]; self.params.repetition()];
        for i in 0..self.params.repetition() {
            for j in 0..batch_count {
                self.oracle.read_monomial_assign(&mut challenge[i][j]);
            }
        }

        for i in 0..self.params.repetition() {
            for j in 0..batch_count {
                for k in 0..self.params.poly_commit_size() {
                    self.params.ringq().mul_add_assign(
                        &challenge[i][j],
                        &open.mask[j][k],
                        &mut open_pf_out.response_mask[i][k],
                    );
                }
                for k in 0..self.params.ajtai_rand_size() {
                    self.params.ringq().mul_add_assign(
                        &challenge[i][j],
                        &open.rand[j][k],
                        &mut open_pf_out.response_rand[i][k],
                    );
                }
            }
        }
    }

    /// Proves the opening without zero-knowledge.
    pub fn prove_opening_plain(&mut self, com: &Commitment, open: &Opening) -> OpeningProof {
        let mut open_pf_out = OpeningProof::new(&self.params);
        self.prove_opening_plain_assign(com, open, &mut open_pf_out);
        return open_pf_out;
    }

    /// Proves the opening without zero-knowledge and writes to open_pf_out.
    pub fn prove_opening_plain_assign(
        &mut self,
        com: &Commitment,
        open: &Opening,
        open_pf_out: &mut OpeningProof,
    ) {
        let batch_count = open.mask.len() - 2;

        self.oracle.reset();
        for i in 0..batch_count {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&com.value[i].value[j]);
            }
        }

        let mut rand_mask = vec![Integer::ZERO; self.params.bigint_commit_size()];
        for i in 0..self.params.repetition() {
            for j in 0..self.params.poly_commit_size() {
                self.uniform_sampler.read_big_mod_assign(&mut rand_mask[j]);
            }
            self.encoder
                .encode_chunk_assign(&rand_mask, &mut open_pf_out.response_mask[i]);
            self.commiter.commit_plain_assign(
                &open_pf_out.response_mask[i],
                &mut open_pf_out.commitment[i],
            );
        }

        for i in 0..self.params.repetition() {
            for j in 0..self.params.ajtai_size() {
                self.oracle.write_poly(&open_pf_out.commitment[i].value[j]);
            }
        }

        let mut challenge =
            vec![vec![self.params.ringq().new_poly(); batch_count]; self.params.repetition()];
        for i in 0..self.params.repetition() {
            for j in 0..batch_count {
                self.oracle.read_monomial_assign(&mut challenge[i][j]);
            }
        }

        for i in 0..self.params.repetition() {
            for j in 0..batch_count {
                for k in 0..self.params.poly_commit_size() {
                    self.params.ringq().mul_add_assign(
                        &challenge[i][j],
                        &open.mask[j][k],
                        &mut open_pf_out.response_mask[i][k],
                    );
                }
            }
        }
    }

    /// Evaluates the polynomial with a proof.
    pub fn evaluate(&mut self, x: Integer, open: &Opening) -> EvaluationProof {
        let mut eval_pf_out = EvaluationProof::new(&self.params);
        self.evaluate_assign(x, open, &mut eval_pf_out);
        return eval_pf_out;
    }

    /// Evaluates the polynomial with a proof and writes to eval_pf_out.
    pub fn evaluate_assign(
        &mut self,
        x: Integer,
        open: &Opening,
        eval_pf_out: &mut EvaluationProof,
    ) {
        let commit_count = open.mask.len() - 2;

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

        for i in 0..self.params.poly_commit_size() {
            eval_pf_out.mask[i].clone_from(&open.mask[commit_count + 1][i]);
            self.params.ringq().mul_add_assign(
                &x_ecd,
                &open.mask[commit_count][i],
                &mut eval_pf_out.mask[i],
            );
            for j in 0..commit_count {
                self.params.ringq().mul_add_assign(
                    &x_pow_ecd[j],
                    &open.mask[j][i],
                    &mut eval_pf_out.mask[i],
                );
            }
        }

        for i in 0..self.params.ajtai_rand_size() {
            eval_pf_out.rand[i].clone_from(&open.rand[commit_count + 1][i]);
            self.params.ringq().mul_add_assign(
                &x_ecd,
                &open.rand[commit_count][i],
                &mut eval_pf_out.rand[i],
            );
            for j in 0..commit_count {
                self.params.ringq().mul_add_assign(
                    &x_pow_ecd[j],
                    &open.rand[j][i],
                    &mut eval_pf_out.rand[i],
                );
            }
        }

        let mut x_pow_basis = vec![Integer::ZERO; self.params.bigint_commit_size()];
        x_pow_buf.assign(1);
        for i in 0..self.params.bigint_commit_size() {
            x_pow_basis[i].assign(&x_pow_buf);
            x_pow_buf *= &x;
            x_pow_buf %= self.params.modulus();
        }

        let mask_dcd = self.encoder.decode_chunk(&eval_pf_out.mask);
        eval_pf_out.value.assign(0);
        let mut prod_buf = Integer::ZERO;
        for i in 0..self.params.bigint_commit_size() {
            prod_buf.assign(&mask_dcd[i]);
            prod_buf *= &x_pow_basis[i];
            eval_pf_out.value += &prod_buf;
        }
        eval_pf_out.value %= self.params.modulus();
    }

    /// Evaluates the polynomial with a proof without hiding.
    pub fn evaluate_plain(&mut self, x: Integer, open: &Opening) -> EvaluationProof {
        let mut eval_pf_out = EvaluationProof::new(&self.params);
        self.evaluate_plain_assign(x, open, &mut eval_pf_out);
        return eval_pf_out;
    }

    /// Evaluates the polynomial with a proof without hiding and writes to eval_pf_out.
    pub fn evaluate_plain_assign(
        &mut self,
        x: Integer,
        open: &Opening,
        eval_pf_out: &mut EvaluationProof,
    ) {
        let commit_count = open.mask.len() - 2;

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

        for i in 0..self.params.poly_commit_size() {
            for j in 0..commit_count {
                self.params.ringq().mul_add_assign(
                    &x_pow_ecd[j],
                    &open.mask[j][i],
                    &mut eval_pf_out.mask[i],
                );
            }
        }

        let mut x_pow_basis = vec![Integer::ZERO; self.params.bigint_commit_size()];
        x_pow_buf.assign(1);
        for i in 0..self.params.bigint_commit_size() {
            x_pow_basis[i].assign(&x_pow_buf);
            x_pow_buf *= &x;
            x_pow_buf %= self.params.modulus();
        }

        let mask_dcd = self.encoder.decode_chunk(&eval_pf_out.mask);
        eval_pf_out.value.assign(0);
        let mut prod_buf = Integer::ZERO;
        for i in 0..self.params.bigint_commit_size() {
            prod_buf.assign(&mask_dcd[i]);
            prod_buf *= &x_pow_basis[i];
            eval_pf_out.value += &prod_buf;
        }
        eval_pf_out.value %= self.params.modulus();
    }
}
