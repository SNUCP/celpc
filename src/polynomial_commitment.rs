use super::csprng::*;
use super::ring::*;
use super::*;
use ethnum::U256;

pub struct PolynomialCommitment {
    pub b0: Vec<Poly>,
    pub b1: Vec<Poly>,
    pub beta0: Vec<Poly>,
    pub beta1: Vec<Poly>,
    pub b0_commit: Vec<Poly>,
    pub b1_commit: Vec<Poly>,

    pub h: Vec<Vec<Poly>>,
    pub eta: Vec<Vec<Poly>>,
    pub h_commit: Vec<Vec<Poly>>,
}

pub struct EvaluationProof {
    pub y_raw: U256,

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
    pub gaussian_sampler: KarneySampler,
    pub oracle: Oracle,

    pub key: &'a CommitKey,

    pub comitter: Comitter<'a>,
}

impl<'a> PolynomialProver<'a> {
    /// Create a new R1CSProver.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> PolynomialProver<'a> {
        PolynomialProver {
            params: params,
            encoder: Encoder::new(params),
            uniform_sampler: UniformSampler::new(),
            gaussian_sampler: KarneySampler::new(),
            oracle: Oracle::new(),

            key: key,

            comitter: Comitter::new(params, key),
        }
    }

    /// Commit creates a polynomial commitment.
    pub fn commit(&mut self, h_raw: &[U256]) -> PolynomialCommitment {
        assert_eq!(h_raw.len(), self.params.m * self.params.n);

        let mut h = vec![vec![self.params.ringq.new_ntt_poly(); self.params.l]; self.params.m];
        for (i, h_raw_chunk) in h_raw.chunks_exact(self.params.n).enumerate() {
            self.encoder
                .encode_randomized_chunk_assign(h_raw_chunk, self.params.s1, &mut h[i]);
        }

        let mut eta = vec![vec![self.params.ringq.new_ntt_poly(); self.params.munu]; self.params.m];
        eta.iter_mut().flatten().for_each(|p| {
            self.gaussian_sampler
                .sample_poly_assign(&self.params.ringq, 0.0, self.params.sig1, p)
        });

        let mut h_commit =
            vec![vec![self.params.ringq.new_ntt_poly(); self.params.mu]; self.params.m];
        for i in 0..self.params.m {
            self.comitter
                .commit_assign(&h[i], &eta[i], &mut h_commit[i]);
        }

        let mut b0_raw = vec![U256::ZERO; self.params.n];
        let mut b1_raw = vec![U256::ZERO; self.params.n];
        for i in 0..self.params.n - 1 {
            let b = self.uniform_sampler.sample_range_u256(self.params.p);
            b0_raw[i + 1] = b.wrapping_neg();
            b1_raw[i] = b;
        }

        let t = ((self.params.m as f64) + 2.0).sqrt();
        let mut b0 = vec![self.params.ringq.new_ntt_poly(); self.params.l];
        self.encoder
            .encode_randomized_chunk_assign(&b0_raw, t * self.params.s3, &mut b0);

        let mut beta0 = vec![self.params.ringq.new_ntt_poly(); self.params.munu];
        beta0.iter_mut().for_each(|p| {
            self.gaussian_sampler.sample_poly_assign(
                &self.params.ringq,
                0.0,
                t * self.params.sig3,
                p,
            )
        });

        let mut b0_commit = vec![self.params.ringq.new_ntt_poly(); self.params.mu];
        self.comitter.commit_assign(&b0, &beta0, &mut b0_commit);

        let mut b1 = vec![self.params.ringq.new_ntt_poly(); self.params.l];
        self.encoder
            .encode_randomized_chunk_assign(&b1_raw, self.params.s1, &mut b1);

        let mut beta1 = vec![self.params.ringq.new_ntt_poly(); self.params.munu];
        beta1.iter_mut().for_each(|p| {
            self.gaussian_sampler
                .sample_poly_assign(&self.params.ringq, 0.0, self.params.sig1, p)
        });

        let mut b1_commit = vec![self.params.ringq.new_ntt_poly(); self.params.mu];
        self.comitter.commit_assign(&b1, &beta1, &mut b1_commit);

        PolynomialCommitment {
            h: h,
            eta: eta,
            h_commit: h_commit,

            b0: b0,
            b1: b1,
            beta0: beta0,
            beta1: beta1,
            b0_commit: b0_commit,
            b1_commit: b1_commit,
        }
    }

    /// Generates an evaluation proof using Polynomial Commitment.
    pub fn evaluate(&mut self, x_raw: U256, pc: PolynomialCommitment) -> EvaluationProof {
        let mut x_buf = self.encoder.encode(&[x_raw]);
        let mut x_raw_powers = vec![U256::ONE; self.params.n];
        for i in 1..self.params.n {
            x_raw_powers[i] = bmod(x_raw_powers[i - 1] * x_raw, self.params.p, self.params.pr);
        }

        let mut e = vec![self.params.ringq.new_ntt_poly(); self.params.l];
        for i in 0..self.params.l {
            self.params.ringq.mul_assign(&x_buf, &pc.b1[i], &mut e[i]);
            self.params.ringq.add_inplace(&pc.b0[i], &mut e[i]);
        }
        let mut eps = vec![self.params.ringq.new_ntt_poly(); self.params.munu];
        for i in 0..self.params.munu {
            self.params
                .ringq
                .mul_assign(&x_buf, &pc.beta1[i], &mut eps[i]);
            self.params.ringq.add_inplace(&pc.beta0[i], &mut eps[i]);
        }

        for i in 0..self.params.m {
            self.encoder
                .encode_assign(&[x_raw_powers[i * self.params.n]], &mut x_buf);
            for j in 0..self.params.l {
                self.params
                    .ringq
                    .mul_add_inplace(&x_buf, &pc.h[i][j], &mut e[j]);
            }
            for j in 0..self.params.munu {
                self.params
                    .ringq
                    .mul_add_inplace(&x_buf, &pc.eta[i][j], &mut eps[j]);
            }
        }

        let mut e_dcd = vec![U256::ZERO; self.params.n];
        self.encoder.decode_chunk_assign(&e, &mut e_dcd);
        let y_raw = inner_product(&x_raw_powers, &e_dcd, self.params.p, self.params.pr);

        EvaluationProof {
            y_raw: y_raw,
            e: e,
            eps: eps,
        }
    }

    pub fn prove(&mut self, pc: PolynomialCommitment) -> OpenProof {
        let k = LAMBDA / (self.params.d.ilog2() as usize);
        let sm = (self.params.m as f64 + 2.0).sqrt();

        let mut g_raw = vec![U256::ZERO; self.params.n];
        let mut g = vec![self.params.ringq.new_ntt_poly(); self.params.l];
        let mut gamma = vec![self.params.ringq.new_ntt_poly(); self.params.munu];
        let mut g_commit = vec![vec![self.params.ringq.new_ntt_poly(); self.params.mu]; k];
        for i in 0..k {
            for j in 0..self.params.n {
                g_raw[j] = self.uniform_sampler.sample_range_u256(self.params.p);
            }
            self.encoder
                .encode_randomized_chunk_assign(&g_raw, self.params.s1, &mut g);
            gamma.iter_mut().for_each(|p| {
                self.gaussian_sampler.sample_poly_assign(
                    &self.params.ringq,
                    0.0,
                    sm * self.params.sig2,
                    p,
                )
            });
            self.comitter.commit_assign(&g, &gamma, &mut g_commit[i]);
        }

        g_commit.iter().flatten().for_each(|p| {
            self.oracle.write_poly(p);
        });
        self.oracle.finalize();

        let mut ch = vec![vec![self.params.ringq.new_ntt_poly(); self.params.m + 1]; k];
        ch.iter_mut()
            .flatten()
            .for_each(|p| self.oracle.read_monomial_assign(&self.params.ringq, p));

        let mut t = vec![vec![self.params.ringq.new_ntt_poly(); self.params.l]; k];
        let mut tau = vec![vec![self.params.ringq.new_ntt_poly(); self.params.munu]; k];
        for i in 0..k {
            for j in 0..self.params.l {
                self.params.ringq.mul_add_inplace(
                    &ch[i][self.params.m + 1],
                    &pc.b1[j],
                    &mut t[i][j],
                );
                self.params.ringq.add_inplace(&g[i], &mut t[i][j]);
            }
            for j in 0..self.params.m {
                for k in 0..self.params.l {
                    self.params
                        .ringq
                        .mul_add_inplace(&ch[i][j], &pc.h[j][k], &mut t[i][k]);
                }
            }

            for j in 0..self.params.munu {
                self.params.ringq.mul_add_inplace(
                    &ch[i][self.params.m + 1],
                    &pc.beta1[j],
                    &mut tau[i][j],
                );
                self.params.ringq.add_inplace(&gamma[j], &mut tau[i][j]);
            }
            for j in 0..self.params.m {
                for k in 0..self.params.munu {
                    self.params
                        .ringq
                        .mul_add_inplace(&ch[i][j], &pc.eta[j][k], &mut tau[i][k]);
                }
            }
        }

        OpenProof {
            ch: ch,
            t: t,
            tau: tau,
        }
    }
}

pub struct PolynomialVerifier<'a> {
    pub params: &'a Parameters,
    pub encoder: Encoder<'a>,
    pub uniform_sampler: UniformSampler,
    pub gaussian_sampler: KarneySampler,
    pub oracle: Oracle,

    pub key: &'a CommitKey,

    pub comitter: Comitter<'a>,

    pub open_bound: f64,
}

impl<'a> PolynomialVerifier<'a> {
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> PolynomialVerifier<'a> {
        PolynomialVerifier {
            params: params,
            encoder: Encoder::new(params),
            uniform_sampler: UniformSampler::new(),
            gaussian_sampler: KarneySampler::new(),
            oracle: Oracle::new(),

            key: key,

            comitter: Comitter::new(params, key),

            open_bound: 0.0,
        }
    }

    pub fn verify_evaluation(
        &mut self,
        x_raw: U256,
        y_raw: U256,
        pc: PolynomialCommitment,
        ep: EvaluationProof,
    ) -> bool {
        let ep_norm2 =
            ep.e.iter()
                .chain(ep.eps.iter())
                .fold(0.0, |acc, p| acc + self.params.ringq.norm(p));
        if ep_norm2 > self.open_bound {
            return false;
        }

        let mut x_raw_powers = vec![U256::ONE; self.params.n];
        for i in 1..self.params.n {
            x_raw_powers[i] = bmod(x_raw_powers[i - 1] * x_raw, self.params.p, self.params.pr);
        }
        let mut e_dcd = vec![U256::ZERO; self.params.n];
        self.encoder.decode_chunk_assign(&ep.e, &mut e_dcd);
        if y_raw != inner_product(&x_raw_powers, &e_dcd, self.params.p, self.params.pr) {
            return false;
        }

        let mut pc_out0 = vec![self.params.ringq.new_ntt_poly(); self.params.l];
        let mut x_buf = self.encoder.encode(&[x_raw]);
        for i in 0..self.params.l {
            self.params
                .ringq
                .mul_add_inplace(&x_buf, &pc.b1[i], &mut pc_out0[i]);
            self.params.ringq.add_inplace(&pc.b0[i], &mut pc_out0[i]);
        }
        for i in 0..self.params.m {
            self.encoder
                .encode_assign(&[x_raw_powers[i * self.params.n]], &mut x_buf);
            for j in 0..self.params.l {
                self.params
                    .ringq
                    .mul_add_inplace(&x_buf, &pc.h[i][j], &mut pc_out0[j]);
            }
        }

        let pc_out1 = self.comitter.commit(&ep.e, &ep.eps);
        if pc_out0.iter().zip(pc_out1.iter()).any(|(a, b)| !a.equal(b)) {
            return false;
        }

        return true;
    }

    pub fn verify(&mut self, pc: PolynomialCommitment, op: OpenProof) -> bool {
        let k = LAMBDA / (self.params.d.ilog2() as usize);

        let mut pc_out = vec![self.params.ringq.new_ntt_poly(); self.params.mu];

        for i in 0..k {
            let op_norm2 = op.t[i]
                .iter()
                .chain(op.tau[i].iter())
                .fold(0.0, |acc, p| acc + self.params.ringq.norm(p));
            if op_norm2 > self.open_bound {
                return false;
            }

            self.comitter
                .commit_assign(&op.t[i], &op.tau[i], &mut pc_out);
            for j in 0..self.params.mu {
                self.params.ringq.mul_sub_inplace(
                    &op.ch[i][self.params.m],
                    &pc.b1_commit[j],
                    &mut pc_out[j],
                );
            }
            for j in 0..self.params.m {
                for k in 0..self.params.mu {
                    self.params.ringq.mul_sub_inplace(
                        &op.ch[i][j],
                        &pc.h_commit[j][k],
                        &mut pc_out[k],
                    );
                }
            }

            pc_out.iter().for_each(|p| self.oracle.write_poly(p));
        }

        let mut m = self.params.ringq.new_ntt_poly();
        self.oracle.finalize();
        if op.ch.iter().flatten().any(|p| {
            self.oracle.read_monomial_assign(&self.params.ringq, &mut m);
            !m.equal(p)
        }) {
            return false;
        }

        return true;
    }
}
