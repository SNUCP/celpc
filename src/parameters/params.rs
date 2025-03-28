use rug::{ops::*, Integer};

use crate::ring::*;

#[derive(Clone)]
pub struct ParametersLiteral {
    /// ajtai_size is the size of the Ajtai Commitment.
    /// Denoted as mu in the paper.
    pub ajtai_size: usize,
    /// ajtai_rand_size is the size of the randomness used in Ajtai Commitment.
    /// Denoted as mu + nu in the paper.
    pub ajtai_rand_size: usize,

    /// degree is the degree of the committing polynomial.
    /// Denoted as N in the paper.
    /// Note that we can commit to a larger/smaller degree polynomial, and this value is just a recommendation.
    pub degree: usize,
    /// bigint_commit_size is the input length of the bigint vector for Ajtai Commitment.
    /// Denoted as n in the paper.
    pub bigint_commit_size: usize,

    /// modulus_base is the base of the modulus of the committing polynomial.
    /// Denoted as b in the paper.
    pub modulus_base: u64,
    /// digits is the exponent of the modulus of the committing polynomial.
    /// Denoted as r in the paper.
    pub digits: usize,

    /// ring_degree is the degree of the internal ring R_q.
    /// Denoted as d in the paper.
    pub ring_degree: usize,
    /// log_ring_modulus is the (log2 value of) modulus of the internal ring R_q.
    /// Denoted as q in the paper.
    pub ring_modulus: Vec<u64>,

    /// commit_std_dev is the standard deviation of randomized encoding used in Ajtai Commitment.
    /// Denoted as s1 in the paper.
    pub commit_std_dev: f64,
    /// opening_proof_std_dev is the standard deviation used in opening proof of Ajtai Commitment.
    /// Denoted as s2 in the paper.
    pub opening_proof_std_dev: f64,
    /// blind_std_dev is the standard deviation of blinding polynomial used in evaluation.
    /// Denoted as s3 in the paper.
    pub blind_std_dev: f64,

    /// commit_rand_std_dev is the standard deviation of error polynomial used in Ajtai Commitment.
    /// Denoted as sigma1 in the paper.
    pub commit_rand_std_dev: f64,
    /// opening_proof_rand_std_dev is the standard deviation used in opening proof of Ajtai Commitment.
    /// Denoted as sigma2 in the paper.
    pub opening_proof_rand_std_dev: f64,
    /// blind_rand_std_dev is the standard deviation of blinding error polynomial used in evaluation.
    /// Denoted as sigma3 in the paper.
    pub blind_rand_std_dev: f64,

    /// open_proof_bound is the bound of opening verification.
    pub open_proof_bound: f64,
    /// eval_bound is the bound of evaluation verification.
    pub eval_bound: f64,
}

#[derive(Clone)]
pub struct Parameters {
    /// ajtai_size is the size of the Ajtai Commitment.
    /// Denoted as mu in the paper.
    ajtai_size: usize,
    /// ajtai_rand_size is the size of the randomness used in Ajtai Commitment.
    /// Denoted as mu + nu in the paper.
    ajtai_rand_size: usize,

    /// degree is the degree of the committing polynomial.
    /// Denoted as N in the paper.
    /// Note that we can commit to a larger/smaller degree polynomial, and this value is just a recommendation.
    degree: usize,
    /// bigint_commit_size is the input length of the bigint vector for Ajtai Commitment.
    /// Denoted as n in the paper.
    bigint_commit_size: usize,
    /// polyVecCommitLength is the input length of the polynomial vector for Ajtai Commitment.
    /// Denoted as l in the paper.
    poly_commit_size: usize,

    // Equals to n / s = n * r / d.
    /// modulus_base is the base of the modulus of the committing polynomial.
    /// Denoted as b in the paper.
    modulus_base: u64,
    /// digits is the exponent of the modulus of the committing polynomial.
    /// Denoted as r in the paper.
    digits: usize,
    /// modulus is the modulus of the committing polynomial.
    /// Denoted as q in the paper.
    modulus: Integer,
    /// slots is the number of "slots" used in high-precision encoding.
    /// Denoted as s in the paper.
    slots: usize,

    /// ringq is the internal ring R_q.
    ringq: Ring,

    /// commit_std_dev is the standard deviation of randomized encoding used in Ajtai Commitment.
    /// Denoted as s1 in the paper.
    commit_std_dev: f64,
    /// opening_proof_std_dev is the standard deviation used in opening proof of Ajtai Commitment.
    /// Denoted as s2 in the paper.
    opening_proof_std_dev: f64,
    /// blind_std_dev is the standard deviation of blinding polynomial used in evaluation.
    /// Denoted as s3 in the paper.
    blind_std_dev: f64,

    /// commit_rand_std_dev is the standard deviation of error polynomial used in Ajtai Commitment.
    /// Denoted as sigma1 in the paper.
    commit_rand_std_dev: f64,
    /// opening_proof_rand_std_dev is the standard deviation used in opening proof of Ajtai Commitment.
    /// Denoted as sigma2 in the paper.
    opening_proof_rand_std_dev: f64,
    /// blind_rand_std_dev is the standard deviation of blinding error polynomial used in evaluation.
    /// Denoted as sigma3 in the paper.
    blind_rand_std_dev: f64,

    /// open_proof_bound is the bound of opening verification.
    open_proof_bound: f64,
    /// eval_bound is the bound of evaluation verification.
    eval_bound: f64,
}

impl Parameters {
    pub fn new(pl: ParametersLiteral) -> Parameters {
        let ringq = Ring::new(pl.ring_degree, &pl.ring_modulus);

        let base_big = Integer::from(pl.modulus_base);
        let modulus = base_big.pow(pl.digits as u32) + 1;

        return Parameters {
            ajtai_size: pl.ajtai_size,
            ajtai_rand_size: pl.ajtai_rand_size,

            degree: pl.degree,
            bigint_commit_size: pl.bigint_commit_size,
            poly_commit_size: pl.bigint_commit_size * pl.digits / pl.ring_degree,

            modulus_base: pl.modulus_base,
            digits: pl.digits,
            modulus: modulus,
            slots: pl.ring_degree / pl.digits,

            ringq: ringq,

            commit_std_dev: pl.commit_std_dev,
            opening_proof_std_dev: pl.opening_proof_std_dev,
            blind_std_dev: pl.blind_std_dev,

            commit_rand_std_dev: pl.commit_rand_std_dev,
            opening_proof_rand_std_dev: pl.opening_proof_rand_std_dev,
            blind_rand_std_dev: pl.blind_rand_std_dev,

            open_proof_bound: pl.open_proof_bound,
            eval_bound: pl.eval_bound,
        };
    }

    /// ajtai_size is the size of the Ajtai Commitment.
    /// Denoted as mu in the paper.
    pub fn ajtai_size(&self) -> usize {
        return self.ajtai_size;
    }

    /// ajtai_rand_size is the size of the randomness used in Ajtai Commitment.
    /// Denoted as mu + nu in the paper.
    pub fn ajtai_rand_size(&self) -> usize {
        self.ajtai_rand_size
    }

    /// degree is the degree of the committing polynomial.
    /// Denoted as N in the paper.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// bigint_commit_size is the input length of the bigint vector for Ajtai Commitment.
    /// Denoted as n in the paper.
    pub fn bigint_commit_size(&self) -> usize {
        self.bigint_commit_size
    }

    /// poly_commit_size is the input length of the polynomial vector for Ajtai Commitment.
    /// Denoted as l in the paper.
    pub fn poly_commit_size(&self) -> usize {
        self.poly_commit_size
    }

    /// modulus_base is the base of the modulus of the committing polynomial.
    /// Denoted as b in the paper.
    pub fn modulus_base(&self) -> u64 {
        self.modulus_base
    }

    /// digits is the exponent of the modulus of the committing polynomial.
    /// Denoted as r in the paper.
    pub fn digits(&self) -> usize {
        self.digits
    }

    /// modulus is the modulus of the committing polynomial.
    /// Denoted as q in the paper.
    pub fn modulus(&self) -> &Integer {
        &self.modulus
    }

    /// slots is the number of "slots" used in high-precision encoding.
    /// Denoted as s in the paper.
    pub fn slots(&self) -> usize {
        self.slots
    }

    /// ringq is the internal ring R_q.
    pub fn ringq(&self) -> &Ring {
        &self.ringq
    }

    /// repetion is the number of repetition count for Proof of Opening Knolwedge protocol.
    pub fn repetition(&self) -> usize {
        return (128. / (1. + (self.ringq.degree() as f64).log2())).ceil() as usize;
    }

    /// commit_std_dev is the standard deviation of randomized encoding used in Ajtai Commitment.
    /// Denoted as s1 in the paper.
    pub fn commit_std_dev(&self) -> f64 {
        self.commit_std_dev
    }

    /// opening_proof_std_dev is the standard deviation used in opening proof of Ajtai Commitment.
    /// Denoted as s2 in the paper.
    pub fn opening_proof_std_dev(&self) -> f64 {
        self.opening_proof_std_dev
    }

    /// blind_std_dev is the standard deviation of blinding polynomial used in evaluation.
    /// Denoted as s3 in the paper.
    pub fn blind_std_dev(&self) -> f64 {
        self.blind_std_dev
    }

    /// commit_rand_std_dev is the standard deviation of error polynomial used in Ajtai Commitment.
    /// Denoted as sigma1 in the paper.
    pub fn commit_rand_std_dev(&self) -> f64 {
        self.commit_rand_std_dev
    }

    /// opening_proof_rand_std_dev is the standard deviation used in opening proof of Ajtai Commitment.
    /// Denoted as sigma2 in the paper.
    pub fn opening_proof_rand_std_dev(&self) -> f64 {
        self.opening_proof_rand_std_dev
    }

    /// blind_rand_std_dev is the standard deviation of blinding error polynomial used in evaluation.
    /// Denoted as sigma3 in the paper.
    pub fn blind_rand_std_dev(&self) -> f64 {
        self.blind_rand_std_dev
    }

    /// open_proof_bound is the bound of opening verification.
    pub fn open_proof_bound(&self) -> f64 {
        self.open_proof_bound
    }

    /// eval_bound is the bound of evaluation verification.
    pub fn eval_bound(&self) -> f64 {
        self.eval_bound
    }
}
