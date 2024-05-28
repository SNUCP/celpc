use super::csprng::*;
use super::ring::*;
use super::*;
use ethnum::U256;
use primitive_types::U512;

pub struct PolynomialCommitment {
    pub h: Vec<Vec<Poly>>,
    pub eta: Vec<Vec<Poly>>,
    pub h_commit: Vec<Vec<Poly>>,
}

pub struct EvaluationProof {
    pub e: Vec<Poly>,
    pub eps: Vec<Poly>,
}

pub struct OpenProof {
    pub ch: Vec<Vec<Poly>>,
    pub t: Vec<Vec<Poly>>,
    pub tau: Vec<Vec<Poly>>,
}

pub struct PolynomialProver<'a> {
    pub params: &'a Parameters,
    pub encoder: Encoder<'a>,
    pub uniform_sampler: UniformSampler,

    pub s1_encoder: EncoderRandSmall<'a>,
    pub s2_encoder: EncoderRandSmall<'a>,
    pub s3_encoder: EncoderRandLarge<'a>,

    pub sig1_sampler: CDTSampler,
    pub sig2_sampler: CDTSampler,
    pub sig3_sampler: ConvolveSampler,

    pub oracle: Oracle,

    pub key: &'a CommitKey,

    pub committer: Committer<'a>,
}

impl<'a> PolynomialProver<'a> {
    /// Create a new R1CSProver.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> PolynomialProver<'a> {
        PolynomialProver {
            params: params,
            encoder: Encoder::new(params),
            uniform_sampler: UniformSampler::new(),
            oracle: Oracle::new(),

            s1_encoder: EncoderRandSmall::new(params, params.s1),
            s2_encoder: EncoderRandSmall::new(params, ((params.m + 2) as f64).sqrt() * params.s2),
            s3_encoder: EncoderRandLarge::new(params, ((params.m + 2) as f64).sqrt() * params.s3),

            sig1_sampler: CDTSampler::new(0.0, params.sig1),
            sig2_sampler: CDTSampler::new(0.0, ((params.m + 2) as f64).sqrt() * params.sig2),
            sig3_sampler: ConvolveSampler::new(((params.m + 2) as f64).sqrt() * params.sig3),

            key: key,

            committer: Committer::new(params, key),
        }
    }

    pub fn commit(&mut self, h_raw: &[U256]) -> PolynomialCommitment {
        assert_eq!(h_raw.len(), self.params.m * self.params.n);

        let mut b0 = vec![U256::ZERO; self.params.n];
        let mut b1 = vec![U256::ZERO; self.params.n];
        for i in 0..self.params.n - 1 {
            b0[i] = self.uniform_sampler.sample_range_u256(self.params.p);
            b1[i + 1] = self.params.p - b0[i];
        }

        let mut h = vec![vec![self.params.ringq.new_poly(); self.params.l]; self.params.m + 2];
        for (i, h_raw_chunk) in h_raw.chunks_exact(self.params.n).enumerate() {
            self.s1_encoder
                .encode_randomized_chunk_assign(h_raw_chunk, &mut h[i]);
        }
        self.s1_encoder
            .encode_randomized_chunk_assign(&b0, &mut h[self.params.m]);
        self.s3_encoder
            .encode_randomized_chunk_assign(&b1, &mut h[self.params.m + 1]);

        let mut eta = vec![vec![self.params.ringq.new_poly(); self.params.munu]; self.params.m + 2];
        for i in 0..self.params.m + 1 {
            for j in 0..self.params.munu {
                self.sig1_sampler
                    .sample_poly_assign(&self.params.ringq, &mut eta[i][j]);
            }
        }
        for j in 0..self.params.munu {
            self.sig3_sampler
                .sample_poly_assign(&self.params.ringq, &mut eta[self.params.m + 1][j]);
        }
        let mut h_commit =
            vec![vec![self.params.ringq.new_ntt_poly(); self.params.mu]; self.params.m + 2];
        for i in 0..self.params.m + 2 {
            self.committer
                .commit_assign(&h[i], &eta[i], &mut h_commit[i]);
        }

        return PolynomialCommitment {
            h: h,
            eta: eta,
            h_commit: h_commit,
        };
    }

    pub fn commit_nozk(&mut self, h_raw: &[U256]) -> PolynomialCommitment {
        assert_eq!(h_raw.len(), self.params.m * self.params.n);

        let mut b0 = vec![U256::ZERO; self.params.n];
        let mut b1 = vec![U256::ZERO; self.params.n];
        for i in 0..self.params.n - 1 {
            b0[i] = self.uniform_sampler.sample_range_u256(self.params.p);
            b1[i + 1] = self.params.p - b0[i];
        }

        let mut h = vec![vec![self.params.ringq.new_poly(); self.params.l]; self.params.m + 2];
        for (i, h_raw_chunk) in h_raw.chunks_exact(self.params.n).enumerate() {
            self.encoder.encode_chunk_assign(h_raw_chunk, &mut h[i]);
        }
        self.encoder.encode_chunk_assign(&b0, &mut h[self.params.m]);
        self.encoder
            .encode_chunk_assign(&b1, &mut h[self.params.m + 1]);

        let eta = vec![vec![self.params.ringq.new_poly(); self.params.munu]; self.params.m + 2];

        let mut h_commit =
            vec![vec![self.params.ringq.new_ntt_poly(); self.params.mu]; self.params.m + 2];
        for i in 0..self.params.m + 2 {
            self.committer.commit_nozk_assign(&h[i], &mut h_commit[i]);
        }

        return PolynomialCommitment {
            h: h,
            eta: eta,
            h_commit: h_commit,
        };
    }

    pub fn evaluate(&mut self, x_raw: U256, pc: &PolynomialCommitment) -> (U256, EvaluationProof) {
        let x = self.encoder.encode(&[x_raw]);

        let x512n_raw = mod_up(mod_exp(x_raw, self.params.n, self.params.p512));
        let mut xn_powers = vec![self.params.ringq.new_ntt_poly(); self.params.m];
        let mut xn_power_raw = U512::one();
        self.encoder
            .encode_assign(&[mod_down(xn_power_raw)], &mut xn_powers[0]);
        for i in 1..self.params.m {
            xn_power_raw = (xn_power_raw * x512n_raw) % self.params.p512;
            self.encoder
                .encode_assign(&[mod_down(xn_power_raw)], &mut xn_powers[i]);
        }

        let mut e = vec![self.params.ringq.new_ntt_poly(); self.params.l];
        for i in 0..self.params.l {
            e[i].coeffs.clone_from(&pc.h[self.params.m + 1][i].coeffs);
            self.params
                .ringq
                .mul_add_inplace(&x, &pc.h[self.params.m][i], &mut e[i]);
            for j in 0..self.params.m {
                self.params
                    .ringq
                    .mul_add_inplace(&xn_powers[j], &pc.h[j][i], &mut e[i]);
            }
        }

        let mut eps = vec![self.params.ringq.new_ntt_poly(); self.params.munu];
        for i in 0..self.params.munu {
            eps[i]
                .coeffs
                .clone_from(&pc.eta[self.params.m + 1][i].coeffs);
            self.params
                .ringq
                .mul_add_inplace(&x, &pc.eta[self.params.m][i], &mut eps[i]);
            for j in 0..self.params.m {
                self.params
                    .ringq
                    .mul_add_inplace(&xn_powers[j], &pc.eta[j][i], &mut eps[i]);
            }
        }

        let mut e_intt = e.clone();
        for p in e_intt.iter_mut() {
            self.params.ringq.intt(p);
        }

        let mut e_dcd = vec![U256::ZERO; self.params.n];
        self.encoder.decode_chunk_assign(&e_intt, &mut e_dcd);

        let mut y_raw = mod_up(e_dcd[self.params.n - 1]);
        let x512_raw = mod_up(x_raw);
        for i in (0..self.params.n - 1).rev() {
            y_raw = y_raw * x512_raw + mod_up(e_dcd[i]);
            y_raw %= self.params.p512;
        }

        return (mod_down(y_raw), EvaluationProof { e: e, eps: eps });
    }

    pub fn prove(&mut self, pc: &PolynomialCommitment) -> OpenProof {
        let kp = ((LAMBDA as f64) / (1.0 + (self.params.d.ilog2() as f64))).ceil() as usize;

        // let s2_m = ((self.params.m + 2) as f64).sqrt() * self.params.s2;
        // let sig2_m = ((self.params.m + 2) as f64).sqrt() * self.params.sig2;

        let mut g_raw = vec![vec![U256::ZERO; self.params.n]; kp];
        let mut g = vec![vec![self.params.ringq.new_ntt_poly(); self.params.l]; kp];
        let mut gamma = vec![vec![self.params.ringq.new_ntt_poly(); self.params.munu]; kp];
        let mut g_commit = vec![vec![self.params.ringq.new_ntt_poly(); self.params.mu]; kp];
        for i in 0..kp {
            g_raw[i].fill_with(|| self.uniform_sampler.sample_range_u256(self.params.p));
            // self.encoder
            //     .encode_randomized_chunk_assign(&g_raw[i], s2_m, &mut g[i]);
            self.s2_encoder
                .encode_randomized_chunk_assign(&g_raw[i], &mut g[i]);
            gamma[i].iter_mut().for_each(|p| {
                // self.gaussian_sampler
                //     .sample_poly_exact_assign(&self.params.ringq, 0, sig2_m, p);
                self.sig2_sampler.sample_poly_assign(&self.params.ringq, p);
            });
            self.committer
                .commit_assign(&g[i], &gamma[i], &mut g_commit[i]);
            g_commit[i].iter().for_each(|p| self.oracle.write_poly(p));
        }

        self.oracle.finalize();
        let mut ch = vec![vec![self.params.ringq.new_ntt_poly(); self.params.m + 1]; kp];
        ch.iter_mut()
            .flatten()
            .for_each(|p| self.oracle.read_monomial_assign(&self.params.ringq, p));

        let mut t = vec![vec![self.params.ringq.new_ntt_poly(); self.params.l]; kp];
        let mut tau = vec![vec![self.params.ringq.new_ntt_poly(); self.params.munu]; kp];
        for i in 0..kp {
            t[i].clone_from(&g[i]);
            for j in 0..self.params.m + 1 {
                for k in 0..self.params.l {
                    self.params
                        .ringq
                        .mul_add_inplace(&ch[i][j], &pc.h[j][k], &mut t[i][k]);
                }
            }

            tau[i].clone_from(&gamma[i]);
            for j in 0..self.params.m + 1 {
                for k in 0..self.params.munu {
                    self.params
                        .ringq
                        .mul_add_inplace(&ch[i][j], &pc.eta[j][k], &mut tau[i][k]);
                }
            }
        }

        return OpenProof {
            ch: ch,
            t: t,
            tau: tau,
        };
    }

    pub fn prove_nozk(&mut self, pc: &PolynomialCommitment) -> OpenProof {
        let kp = ((LAMBDA as f64) / (1.0 + (self.params.d.ilog2() as f64))).ceil() as usize;

        let mut g_raw = vec![vec![U256::ZERO; self.params.n]; kp];
        let mut g = vec![vec![self.params.ringq.new_ntt_poly(); self.params.l]; kp];
        let mut g_commit = vec![vec![self.params.ringq.new_ntt_poly(); self.params.mu]; kp];
        for i in 0..kp {
            g_raw[i].fill_with(|| self.uniform_sampler.sample_range_u256(self.params.p));
            self.encoder.encode_chunk_assign(&g_raw[i], &mut g[i]);
            self.committer.commit_nozk_assign(&g[i], &mut g_commit[i]);
            g_commit[i].iter().for_each(|p| self.oracle.write_poly(p));
        }

        self.oracle.finalize();
        let mut ch = vec![vec![self.params.ringq.new_ntt_poly(); self.params.m + 1]; kp];
        ch.iter_mut()
            .flatten()
            .for_each(|p| self.oracle.read_monomial_assign(&self.params.ringq, p));

        let mut t = vec![vec![self.params.ringq.new_ntt_poly(); self.params.l]; kp];
        let tau = vec![vec![self.params.ringq.new_ntt_poly(); self.params.munu]; kp];
        for i in 0..kp {
            t[i].clone_from(&g[i]);
            for j in 0..self.params.m + 1 {
                for k in 0..self.params.l {
                    self.params
                        .ringq
                        .mul_add_inplace(&ch[i][j], &pc.h[j][k], &mut t[i][k]);
                }
            }
        }

        return OpenProof {
            ch: ch,
            t: t,
            tau: tau,
        };
    }
}

pub struct PolynomialVerifier<'a> {
    pub params: &'a Parameters,
    pub encoder: Encoder<'a>,
    pub uniform_sampler: UniformSampler,
    pub oracle: Oracle,

    pub key: &'a CommitKey,

    pub committer: Committer<'a>,
}

impl<'a> PolynomialVerifier<'a> {
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> PolynomialVerifier<'a> {
        PolynomialVerifier {
            params: params,
            encoder: Encoder::new(params),
            uniform_sampler: UniformSampler::new(),
            oracle: Oracle::new(),

            key: key,

            committer: Committer::new(params, key),
        }
    }

    pub fn verify_evaluation(
        &mut self,
        x_raw: U256,
        y_raw: U256,
        pc: &PolynomialCommitment,
        ep: &EvaluationProof,
    ) -> bool {
        let mut e_intt = ep.e.clone();
        for p in e_intt.iter_mut() {
            self.params.ringq.intt(p);
        }
        let mut eps_intt = ep.eps.clone();
        for p in eps_intt.iter_mut() {
            self.params.ringq.intt(p);
        }

        let ev_norm2 = e_intt
            .iter()
            .chain(eps_intt.iter())
            .fold(U256::ZERO, |acc, p| acc + self.params.ringq.norm(p));

        if U256log2(ev_norm2) > 2.0 * self.params.log_bound_eval {
            return false;
        }

        let x = self.encoder.encode(&[x_raw]);
        let x512n_raw = mod_up(mod_exp(x_raw, self.params.n, self.params.p512));
        let mut xn_powers = vec![self.params.ringq.new_ntt_poly(); self.params.m];
        let mut xn_power_raw = U512::one();
        self.encoder
            .encode_assign(&[mod_down(xn_power_raw)], &mut xn_powers[0]);
        for i in 1..self.params.m {
            xn_power_raw = (xn_power_raw * x512n_raw) % self.params.p512;
            self.encoder
                .encode_assign(&[mod_down(xn_power_raw)], &mut xn_powers[i]);
        }

        let mut e_dcd = vec![U256::ZERO; self.params.n];
        self.encoder.decode_chunk_assign(&e_intt, &mut e_dcd);

        let mut y_raw_check = mod_up(e_dcd[self.params.n - 1]);
        let x512_raw = mod_up(x_raw);
        for i in (0..self.params.n - 1).rev() {
            y_raw_check = y_raw_check * x512_raw + mod_up(e_dcd[i]);
            y_raw_check %= self.params.p512;
        }
        if y_raw != mod_down(y_raw_check) {
            return false;
        }

        let mut commit_check = self.committer.commit(&ep.e, &ep.eps);
        for i in 0..self.params.mu {
            self.params
                .ringq
                .sub_inplace(&pc.h_commit[self.params.m + 1][i], &mut commit_check[i]);
            self.params.ringq.mul_sub_inplace(
                &x,
                &pc.h_commit[self.params.m][i],
                &mut commit_check[i],
            );

            for j in 0..self.params.m {
                self.params.ringq.mul_sub_inplace(
                    &xn_powers[j],
                    &pc.h_commit[j][i],
                    &mut commit_check[i],
                );
            }
        }

        if commit_check.iter().any(|p| !p.is_zero()) {
            return false;
        }

        return true;
    }

    pub fn verify(&mut self, pc: &PolynomialCommitment, op: &OpenProof) -> bool {
        let kp = ((LAMBDA as f64) / (1.0 + (self.params.d.ilog2() as f64))).ceil() as usize;

        for i in 0..kp {
            let op_norm2 = op.t[i]
                .iter()
                .chain(op.tau[i].iter())
                .fold(U256::ZERO, |acc, p| acc + self.params.ringq.norm(p));

            if U256log2(op_norm2) > 2.0 * self.params.log_bound_open {
                return false;
            }

            let mut g_commit_check = vec![self.params.ringq.new_ntt_poly(); self.params.mu];
            self.committer
                .commit_assign(&op.t[i], &op.tau[i], &mut g_commit_check);
            for j in 0..self.params.m + 1 {
                for k in 0..self.params.mu {
                    self.params.ringq.mul_sub_inplace(
                        &op.ch[i][j],
                        &pc.h_commit[j][k],
                        &mut g_commit_check[k],
                    );
                }
            }
            g_commit_check
                .iter()
                .for_each(|p| self.oracle.write_poly(p));
        }

        let mut mono_check = self.params.ringq.new_ntt_poly();
        self.oracle.finalize();
        for i in 0..kp {
            for j in 0..self.params.m + 1 {
                self.oracle
                    .read_monomial_assign(&self.params.ringq, &mut mono_check);
                if !op.ch[i][j].equal(&mono_check) {
                    return false;
                }
            }
        }

        return true;
    }
}
