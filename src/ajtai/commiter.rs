use super::entities::*;
use crate::{parameters::Parameters, ring::*};

/// AjtaiCommiter commits an Ajtai Commitment.
#[derive(Clone)]
pub struct AjtaiCommiter {
    parameters: Parameters,
    commit_key: AjtaiCommitKey,
}

impl AjtaiCommiter {
    /// Creates a new AjtaiCommiter.
    pub fn new(params: &Parameters, commit_key: &AjtaiCommitKey) -> AjtaiCommiter {
        AjtaiCommiter {
            parameters: params.clone(),
            commit_key: commit_key.clone(),
        }
    }

    /// Commits p using opening r.
    pub fn commit(&self, p: &[Poly], r: &[Poly]) -> AjtaiCommitment {
        let mut com_out = AjtaiCommitment::new(&self.parameters);
        self.commit_assign(p, r, &mut com_out);
        com_out
    }

    /// Commits p using opening r and writes it to com_out.
    pub fn commit_assign(&self, p: &[Poly], r: &[Poly], com_out: &mut AjtaiCommitment) {
        for i in 0..self.parameters.ajtai_size() {
            self.parameters.ringq().mul_assign(
                &self.commit_key.a0[i][0],
                &p[0],
                &mut com_out.value[i],
            );
            for j in 1..self.parameters.poly_commit_size() {
                self.parameters.ringq().mul_add_assign(
                    &self.commit_key.a0[i][j],
                    &p[j],
                    &mut com_out.value[i],
                );
            }

            for j in 0..self.parameters.ajtai_rand_size() - self.parameters.ajtai_size() {
                self.parameters.ringq().mul_add_assign(
                    &self.commit_key.a1[i][j],
                    &r[j],
                    &mut com_out.value[i],
                );
            }
            self.parameters.ringq().add_inplace(
                &r[self.parameters.ajtai_rand_size() - self.parameters.ajtai_size() + i],
                &mut com_out.value[i],
            );
        }
    }

    /// Commits p without hiding.
    pub fn commit_plain(&self, p: &[Poly]) -> AjtaiCommitment {
        let mut com_out = AjtaiCommitment::new(&self.parameters);
        self.commit_plain_assign(p, &mut com_out);
        com_out
    }

    /// Commits p without hiding, and writes it to com_out.
    pub fn commit_plain_assign(&self, p: &[Poly], com_out: &mut AjtaiCommitment) {
        for i in 0..self.parameters.ajtai_size() {
            self.parameters.ringq().mul_assign(
                &self.commit_key.a0[i][0],
                &p[0],
                &mut com_out.value[i],
            );
            for j in 1..self.parameters.poly_commit_size() {
                self.parameters.ringq().mul_add_assign(
                    &self.commit_key.a0[i][j],
                    &p[j],
                    &mut com_out.value[i],
                );
            }
        }
    }
}
