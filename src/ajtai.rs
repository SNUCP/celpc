use crate::{csprng::*, ring::*, *};

pub struct AjtaiProof {
    // Length [4k+2][mu]
    pub ncommit: Vec<Vec<Poly>>,
    // Length [4k+2][l]
    pub m: Vec<Vec<Poly>>,
    // Length [4k+2][munu]
    pub mu: Vec<Vec<Poly>>,
}

pub struct AjtaiProver<'a> {
    pub params: &'a Parameters,
    pub sampler: KarneySampler,
    pub oracle: Oracle,
    pub comitter: Comitter<'a>,
    pub key: &'a CommitKey,

    pub msg_stddev: f64,
    pub rand_stddev: f64,
}

impl<'a> AjtaiProver<'a> {
    /// Create a new AjtaiProver.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> AjtaiProver<'a> {
        let stddev = (4.0 * (params.k as f64) + 2.0).sqrt();
        AjtaiProver {
            params: params,
            sampler: KarneySampler::new(),
            oracle: Oracle::new(),
            comitter: Comitter::new(params, key),
            key: key,

            msg_stddev: stddev * params.s1,
            rand_stddev: stddev * params.sig1,
        }
    }

    pub fn prove(
        &mut self,
        m: &[Vec<Poly>],
        mu: &[Vec<Poly>],
        mcommit: &[Vec<Poly>],
    ) -> AjtaiProof {
        let params = self.params;

        // Step 1
        let mut n = vec![vec![params.ringq.new_ntt_poly(); params.l]; params.r];
        let mut nu = vec![vec![params.ringq.new_ntt_poly(); params.munu]; params.r];
        let mut ncommit = vec![vec![params.ringq.new_ntt_poly(); params.mu]; params.r];
        for i in 0..params.r {
            for j in 0..params.l {
                self.sampler
                    .sample_poly_assign(&params.ringq, 0.0, self.msg_stddev, &mut n[i][j]);
            }

            for j in 0..params.munu {
                self.sampler.sample_poly_assign(
                    &params.ringq,
                    0.0,
                    self.rand_stddev,
                    &mut nu[i][j],
                );
            }

            self.comitter.commit_assign(&n[i], &nu[i], &mut ncommit[i]);
        }

        // Step 2
        for i in 0..4 * params.k + 2 {
            for j in 0..params.mu {
                self.oracle.write_poly(&mcommit[i][j]);
            }
        }
        for i in 0..params.r {
            for j in 0..params.mu {
                self.oracle.write_poly(&ncommit[i][j]);
            }
        }
        self.oracle.finalize();

        let mut ch = vec![vec![params.ringq.new_ntt_poly(); 4 * params.k + 2]; params.r];
        for i in 0..params.r {
            for j in 0..4 * params.k + 2 {
                self.oracle
                    .read_monomial_assign(&params.ringq, &mut ch[i][j]);
            }
        }

        // Step 3
        let mut mout = vec![vec![params.ringq.new_ntt_poly(); params.l]; params.r];
        let mut muout = vec![vec![params.ringq.new_ntt_poly(); params.munu]; params.r];
        for i in 0..params.r {
            for j in 0..params.l {
                mout[i][j].clone_from(&n[i][j]);
                for k in 0..4 * params.k + 2 {
                    params
                        .ringq
                        .mul_add_inplace(&ch[i][k], &m[k][j], &mut mout[i][j]);
                }
            }

            for j in 0..params.munu {
                muout[i][j].clone_from(&nu[i][j]);
                for k in 0..4 * params.k + 2 {
                    params
                        .ringq
                        .mul_add_inplace(&ch[i][k], &mu[k][j], &mut muout[i][j]);
                }
            }
        }

        AjtaiProof {
            ncommit: ncommit,
            m: mout,
            mu: muout,
        }
    }
}

pub struct AjtaiVerifier<'a> {
    pub params: &'a Parameters,
    pub oracle: Oracle,
    pub comitter: Comitter<'a>,
    pub key: &'a CommitKey,

    mbound2: f64,
    mboundinf: f64,
    mubound2: f64,
    muboundinf: f64,
}

impl<'a> AjtaiVerifier<'a> {
    /// Create a new AjtaiVerifier.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> AjtaiVerifier<'a> {
        let ak = 4.0 * (params.k as f64) + 2.0;
        let bf = params.b as f64;

        AjtaiVerifier {
            params: params,
            oracle: Oracle::new(),
            comitter: Comitter::new(params, key),
            key: key,

            mbound2: ((bf + 1.0) * params.s1 + params.s2)
                * (ak * (params.l * params.n) as f64).sqrt(),
            mboundinf: 5.0 * ((bf + 1.0) * params.s1 + params.s2) * ak.sqrt(),
            mubound2: (params.sig1 + params.sig2) * (ak * (params.munu * params.n) as f64).sqrt(),
            muboundinf: 5.0 * (params.sig1 + params.sig2) * ak.sqrt(),
        }
    }

    pub fn verify(&mut self, mcommit: &[Vec<Poly>], proof: &AjtaiProof) -> bool {
        let params = self.params;

        for i in 0..4 * params.k + 2 {
            for j in 0..params.mu {
                self.oracle.write_poly(&mcommit[i][j]);
            }
        }
        for i in 0..params.r {
            for j in 0..params.mu {
                self.oracle.write_poly(&proof.ncommit[i][j]);
            }
        }
        self.oracle.finalize();

        let mut ch = vec![vec![params.ringq.new_ntt_poly(); 4 * params.k + 2]; params.r];
        for i in 0..params.r {
            for j in 0..4 * params.k + 2 {
                self.oracle
                    .read_monomial_assign(&params.ringq, &mut ch[i][j]);
            }
        }

        for i in 0..params.r {
            let (l2, linf) = proof.m[i].iter().fold((0.0, 0.0), |acc, p| {
                let (l2, linf) = params.ringq.norm(p);
                (acc.0 + l2, f64::max(acc.1, linf))
            });
            if !(l2.sqrt() < self.mbound2 && linf < self.mboundinf) {
                return false;
            }

            let (l2, linf) = proof.mu[i].iter().fold((0.0, 0.0), |acc, p| {
                let (l2, linf) = params.ringq.norm(p);
                (acc.0 + l2, f64::max(acc.1, linf))
            });
            if !(l2.sqrt() < self.mubound2 && linf < self.muboundinf) {
                return false;
            }

            let commit_lhs = self.comitter.commit(&proof.m[i], &proof.mu[i]);
            let mut commit_rhs = params.ringq.new_ntt_poly();
            for j in 0..params.mu {
                commit_rhs.clone_from(&proof.ncommit[i][j]);
                for k in 0..4 * params.k + 2 {
                    params
                        .ringq
                        .mul_add_inplace(&ch[i][k], &mcommit[k][j], &mut commit_rhs);
                }

                if !commit_lhs[j].equal(&commit_rhs) {
                    return false;
                }
            }
        }

        return true;
    }
}
