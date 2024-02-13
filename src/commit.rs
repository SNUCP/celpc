use crate::{csprng::*, ring::*, Parameters};

#[derive(Debug, Clone)]
pub struct CommitKey {
    /// A0 is of size [mu][l].
    pub A0: Vec<Vec<Poly>>,
    /// A1 is of size [mu][nu].
    pub A1: Vec<Vec<Poly>>,
}

impl CommitKey {
    /// Create a new commit key.
    pub fn new(params: &Parameters, seed: &[u8]) -> CommitKey {
        let mut s = UniformSampler::new_with_seed(seed);

        let mut A0 = vec![vec![params.ringq.new_ntt_poly(); params.l]; params.mu];
        let mut A1 = vec![vec![params.ringq.new_ntt_poly(); params.nu]; params.mu];

        for i in 0..params.mu {
            for j in 0..params.l {
                s.sample_poly_assign(&params.ringq, &mut A0[i][j]);
            }
            for j in 0..params.nu {
                s.sample_poly_assign(&params.ringq, &mut A1[i][j]);
            }
        }

        CommitKey { A0, A1 }
    }
}

pub struct Committer<'a> {
    pub params: &'a Parameters,
    pub key: &'a CommitKey,
}

impl<'a> Committer<'a> {
    /// Create a new commiter.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> Committer<'a> {
        Committer {
            params: params,
            key: key,
        }
    }

    /// Commits a new packed message with randomness rand.
    /// Message should be of size l, which is a randomized encoding of n Z_p elements.
    /// Randomness should be of size munu.
    /// Commitment is of size mu.
    pub fn commit(&self, message: &[Poly], rand: &[Poly]) -> Vec<Poly> {
        let mut commit_out = vec![self.params.ringq.new_ntt_poly(); self.params.mu];
        self.commit_assign(message, rand, &mut commit_out);
        commit_out
    }

    /// Commits a new packed message with randomness r.
    /// Message should be of size l, which is a randomized encoding of n Z_p elements.
    /// Randomness should be of size munu.
    /// Commitment should be of size mu.
    pub fn commit_assign(self: &Self, message: &[Poly], rand: &[Poly], commit_out: &mut [Poly]) {
        let params = self.params;

        for i in 0..params.mu {
            commit_out[i].is_ntt = true;
            commit_out[i].clear();
        }

        // A0 * m
        for i in 0..params.mu {
            for j in 0..params.l {
                params
                    .ringq
                    .mul_add_inplace(&self.key.A0[i][j], &message[j], &mut commit_out[i]);
            }
        }

        // A1 * mu
        for i in 0..params.mu {
            for j in 0..params.nu {
                params
                    .ringq
                    .mul_add_inplace(&self.key.A1[i][j], &rand[j], &mut commit_out[i]);
            }
            params
                .ringq
                .add_inplace(&rand[params.mu + i], &mut commit_out[i]);
        }
    }

    /// Commits a new packed message with randomness r.
    /// Message should be of size l, which is a randomized encoding of n Z_p elements.
    /// Commitment should be of size mu.
    pub fn commit_nozk_assign(self: &Self, message: &[Poly], commit_out: &mut [Poly]) {
        let params = self.params;

        for i in 0..params.mu {
            commit_out[i].is_ntt = true;
            commit_out[i].clear();
        }

        // A0 * m
        for i in 0..params.mu {
            for j in 0..params.l {
                params
                    .ringq
                    .mul_add_inplace(&self.key.A0[i][j], &message[j], &mut commit_out[i]);
            }
        }
    }
}
