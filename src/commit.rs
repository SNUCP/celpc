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

pub struct Comitter<'a> {
    pub params: &'a Parameters,
    pub key: &'a CommitKey,
}

impl<'a> Comitter<'a> {
    /// Create a new commiter.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> Comitter<'a> {
        Comitter {
            params: params,
            key: key,
        }
    }

    /// Commits a new packed message with randomness r.
    /// Message is of size l.
    /// Randomness is of size munu.
    /// Commitment is of size mu.
    pub fn commit(&self, message: &[Poly], mu: &[Poly]) -> Vec<Poly> {
        let mut commit_out = vec![self.params.ringq.new_ntt_poly(); self.params.mu];
        self.commit_assign(message, mu, &mut commit_out);
        commit_out
    }

    /// Commits a new packed message with randomness r.
    /// Message is of size l.
    /// Randomness is of size munu.
    /// Commitment is of size mu.
    pub fn commit_assign(self: &Self, m: &[Poly], mu: &[Poly], commit_out: &mut [Poly]) {
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
                    .mul_add_inplace(&self.key.A0[i][j], &m[j], &mut commit_out[i]);
            }
        }

        // A1 * mu
        for i in 0..params.mu {
            for j in 0..params.nu {
                params
                    .ringq
                    .mul_add_inplace(&self.key.A1[i][j], &mu[j], &mut commit_out[i]);
            }
            params
                .ringq
                .add_inplace(&mu[params.mu + i], &mut commit_out[i]);
        }
    }

    /// Opens a commitment.
    /// Message is of size l.
    /// Randomness is of size munu.
    /// Commitment is of size mu.
    pub fn open(&self, m: &[Poly], mu: &[Poly], bound: (f64, f64), commit: &[Poly]) -> bool {
        let params = &self.params;

        let (l2, linf) = m.iter().chain(mu).fold((0.0, 0.0), |acc, x| {
            let (l2, linf) = params.ringq.norm(x);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < bound.0 && linf < bound.1) {
            return false;
        }

        let commit_out = self.commit(m, mu);
        for i in 0..params.mu {
            if !&commit[i].equal(&commit_out[i]) {
                return false;
            }
        }

        return true;
    }
}

#[cfg(test)]
mod tests {
    use crate::csprng::*;
    use crate::*;
    use primitive_types::U256;

    #[test]
    fn test_commit() {
        let params = Parameters::small();
        let ecd = Encoder::new(&params);
        let mut us = UniformSampler::new();
        let mut gs = KarneySampler::new();

        let mut m = vec![params.ringq.new_ntt_poly(); params.l];
        let mut msg = vec![U256::from(0); params.m];
        for i in 0..params.l {
            for j in 0..params.m {
                msg[j] = us.sample_u256() % params.p;
            }
            ecd.encode_assign(&msg, &mut m[i]);
        }

        let mut mu = vec![params.ringq.new_ntt_poly(); params.munu];
        for i in 0..params.munu {
            gs.sample_poly_assign(&params.ringq, 0.0, 3.2, &mut mu[i]);
        }

        let ck = CommitKey::new(&params, &[0u8; 32]);
        let commiter = Comitter::new(&params, &ck);
        let commit = commiter.commit(&m, &mu);
        if !commiter.open(&m, &mu, (f64::exp2(30.0), f64::exp(30.0)), &commit) {
            panic!("commitment failed to open");
        }
    }
}
