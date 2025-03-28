use rug::Integer;

use crate::{ajtai::AjtaiCommitment, parameters::Parameters, ring::Poly};

/// Commitment is a polynomial commitment.
pub struct Commitment {
    /// value is a commitment value.
    /// Denoted as h.
    pub value: Vec<AjtaiCommitment>,
}

/// Opening is an opening of a polynomial commitment.
pub struct Opening {
    /// mask is the masking polynomial.
    /// Denoted as h.
    pub mask: Vec<Vec<Poly>>,
    /// rand is the randomness.
    /// Denoted as eta.
    pub rand: Vec<Vec<Poly>>,
}

/// OpeningProof is a proof of an opening.
pub struct OpeningProof {
    /// commitment is the first move of the opening proof.
    pub commitment: Vec<AjtaiCommitment>,
    /// response_mask is the masking in the last move of the opening proof.
    pub response_mask: Vec<Vec<Poly>>,
    /// response_rand is the randomness in the last move of the opening proof.
    pub response_rand: Vec<Vec<Poly>>,
}

/// EvaluationProof is a proof of an evaluation.
pub struct EvaluationProof {
    /// value is the evaluated value of the polynomial.
    pub value: Integer,
    /// mask is the masking polynomial.
    /// Denoted as e.
    pub mask: Vec<Poly>,
    /// rand is the randomness.
    /// Denoted as rho.
    pub rand: Vec<Poly>,
}

impl Commitment {
    /// Creates a new commitment.
    pub fn new(params: &Parameters, degree: usize) -> Commitment {
        let split_count = degree / params.bigint_commit_size();

        return Commitment {
            value: vec![AjtaiCommitment::new(params); split_count + 2],
        };
    }
}

impl Opening {
    /// Creates a new opening.
    pub fn new(params: &Parameters, degree: usize) -> Opening {
        let split_count = degree / params.bigint_commit_size();

        return Opening {
            mask: vec![vec![params.ringq().new_poly(); params.poly_commit_size()]; split_count + 2],
            rand: vec![vec![params.ringq().new_poly(); params.ajtai_rand_size()]; split_count + 2],
        };
    }
}

impl OpeningProof {
    /// Creates a new opening proof.
    pub fn new(params: &Parameters) -> OpeningProof {
        return OpeningProof {
            commitment: vec![AjtaiCommitment::new(params); params.repetition()],
            response_mask: vec![
                vec![params.ringq().new_poly(); params.poly_commit_size()];
                params.repetition()
            ],
            response_rand: vec![
                vec![params.ringq().new_poly(); params.ajtai_rand_size()];
                params.repetition()
            ],
        };
    }
}

impl EvaluationProof {
    /// Creates a new evaluation proof.
    pub fn new(params: &Parameters) -> EvaluationProof {
        return EvaluationProof {
            value: Integer::new(),
            mask: vec![params.ringq().new_poly(); params.poly_commit_size()],
            rand: vec![params.ringq().new_poly(); params.ajtai_rand_size()],
        };
    }
}
