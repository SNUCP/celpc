use crate::{csprng, parameters::Parameters, ring::*};

/// AjtaiCommitKey is a public CRS for Ajtai Commitment.
#[derive(Clone)]
pub struct AjtaiCommitKey {
    /// a0 has dimension ajtai_size * poly_commit_size.
    pub a0: Vec<Vec<Poly>>,
    /// a1 has dimension ajtai_size * (ajtai_rand_size - ajtai_size).
    pub a1: Vec<Vec<Poly>>,
}

impl AjtaiCommitKey {
    /// Creates and Generates a new AjtaiCommitKey.
    pub fn new(params: &Parameters) -> AjtaiCommitKey {
        let mut s = csprng::UniformSampler::new(params);

        let mut a0 =
            vec![vec![params.ringq().new_poly(); params.poly_commit_size()]; params.ajtai_size()];
        for i in 0..params.ajtai_size() {
            for j in 0..params.poly_commit_size() {
                s.read_poly_assign(&mut a0[i][j]);
            }
        }

        let mut a1 =
            vec![
                vec![params.ringq().new_poly(); params.ajtai_rand_size() - params.ajtai_size()];
                params.ajtai_size()
            ];
        for i in 0..params.ajtai_size() {
            for j in 0..params.ajtai_rand_size() - params.ajtai_size() {
                s.read_poly_assign(&mut a1[i][j]);
            }
        }

        AjtaiCommitKey { a0, a1 }
    }
}

/// AjtaiCommitment is a Ajtai AjtaiCommitment.
#[derive(Clone)]
pub struct AjtaiCommitment {
    /// Value has length ajtai_size.
    pub value: Vec<Poly>,
}

impl AjtaiCommitment {
    /// Creates a new Commitment.
    pub fn new(params: &Parameters) -> AjtaiCommitment {
        let com = vec![params.ringq().new_poly(); params.ajtai_size()];

        AjtaiCommitment { value: com }
    }

    /// Checks if two AjtaiCommitments are equal.
    pub fn equals(&self, other: &AjtaiCommitment) -> bool {
        if self.value.len() != other.value.len() {
            return false;
        }

        for i in 0..self.value.len() {
            if !self.value[i].equal(&other.value[i]) {
                return false;
            }
        }

        return true;
    }
}
