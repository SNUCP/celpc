use crate::{csprng::*, parameters::Parameters, ring::*};
use rug::{Assign, Float, Integer};

/// Encoder computes the encoding of a polynomial.
#[derive(Clone)]
pub struct Encoder {
    pub params: Parameters,

    pub gaussian_sampler: COSACSampler,

    pub base_big: Integer,
    pub delta_inv: Vec<f64>,
}

impl Encoder {
    /// Creates a new encoder.
    pub fn new(params: &Parameters) -> Self {
        let base_big = Integer::from(params.modulus_base());

        let prec = params.modulus().significant_bits();
        let base_float: Float = Float::with_val(prec, params.modulus_base());
        let mut modulus_inv: Float = -1 / Float::with_val(prec, params.modulus());
        let mut delta_inv = vec![0.0f64; params.digits()];
        for i in 0..delta_inv.len() {
            delta_inv[i] = modulus_inv.to_f64();
            modulus_inv *= &base_float;
        }

        Self {
            params: params.clone(),
            gaussian_sampler: COSACSampler::new(params),
            base_big,
            delta_inv,
        }
    }

    /// Encodes a bigint vector to polynomial and returns it.
    pub fn encode(&self, v: &[Integer]) -> Poly {
        let mut p_out = self.params.ringq().new_poly();
        self.encode_assign(v, &mut p_out);
        return p_out;
    }

    /// Encodes a bigint vector to polynomial and assigns to p_out.
    pub fn encode_assign(&self, v: &[Integer], p_out: &mut Poly) {
        p_out.clear();

        let mut coeff = Integer::ZERO;
        let mut rem = Integer::ZERO;
        for i in 0..v.len() {
            coeff.assign(&v[i]);
            coeff.modulo_mut(&self.params.modulus());
            for j in 0..(self.params.digits() - 1) {
                rem.assign(&self.base_big);
                coeff.div_rem_mut(&mut rem);
                let r = rem.to_u64().unwrap();
                for k in 0..self.params.ringq().level() {
                    p_out.coeffs[k][j * self.params.slots() + i] = r;
                }
            }
            let r = coeff.to_u64().unwrap();
            for k in 0..self.params.ringq().level() {
                p_out.coeffs[k][(self.params.digits() - 1) * self.params.slots() + i] = r;
            }
        }

        self.params.ringq().ntt_inplace(p_out);
    }

    /// Encodes a bigint vector to multiple polynomials and returns it.
    /// The number of polynomials is v.len() / slots.
    pub fn encode_chunk(&self, v: &[Integer]) -> Vec<Poly> {
        let mut p_out = vec![self.params.ringq().new_poly(); v.len() / self.params.slots()];
        self.encode_chunk_assign(v, &mut p_out);
        return p_out;
    }

    /// Encodes a bigint vector to multiple polynomials and assigns to p_out.
    /// The number of polynomials is v.len() / slots.
    pub fn encode_chunk_assign(&self, v: &[Integer], p_out: &mut [Poly]) {
        for i in 0..v.len() / self.params.slots() {
            self.encode_assign(
                &v[i * self.params.slots()..(i + 1) * self.params.slots()],
                &mut p_out[i],
            );
        }
    }

    /// Encodes a bigint vector to polynomial with randomization and returns it.
    pub fn encode_randomize(&mut self, v: &[Integer], std_dev: f64) -> Poly {
        let mut p_out = self.params.ringq().new_poly();
        self.encode_randomize_assign(v, std_dev, &mut p_out);
        return p_out;
    }

    /// Encodes a bigint vector to polynomial with randomization and assigns to p_out.
    pub fn encode_randomize_assign(&mut self, v: &[Integer], std_dev: f64, p_out: &mut Poly) {
        let mut fp_ecd = vec![0.0f64; self.params.ringq().degree()];

        let mut coeff = Integer::ZERO;
        let mut rem = Integer::ZERO;
        for i in 0..v.len() {
            coeff.assign(&v[i]);
            coeff.modulo_mut(&self.params.modulus());
            for j in 0..(self.params.digits() - 1) {
                rem.assign(&self.base_big);
                coeff.div_rem_mut(&mut rem);
                let r = rem.to_f64();
                fp_ecd[j * self.params.slots() + i] = r;
            }
            let r = coeff.to_f64();
            fp_ecd[(self.params.digits() - 1) * self.params.slots() + i] = r;
        }

        let mut fp_sample = vec![0.0f64; self.params.ringq().degree()];
        for i in 0..self.params.digits() {
            let d = self.params.ringq().degree() - (i + 1) * self.params.slots();
            for j in 0..self.params.ringq().degree() - d {
                fp_sample[j + d] += self.delta_inv[i] * fp_ecd[j];
            }
            for j in self.params.ringq().degree() - d..self.params.ringq().degree() {
                fp_sample[j + d - self.params.ringq().degree()] -= self.delta_inv[i] * fp_ecd[j];
            }
        }

        for i in 0..self.params.ringq().degree() {
            fp_sample[i] = self.gaussian_sampler.read_coset_f64(fp_sample[i], std_dev);
        }

        let base_f64 = self.params.modulus_base() as f64;
        for i in 0..self.params.ringq().degree() - self.params.slots() {
            fp_ecd[i + self.params.slots()] =
                fp_sample[i] - base_f64 * fp_sample[i + self.params.slots()];
        }
        for i in self.params.ringq().degree() - self.params.slots()..self.params.ringq().degree() {
            fp_ecd[i + self.params.slots() - self.params.ringq().degree()] = -fp_sample[i]
                - base_f64 * fp_sample[i + self.params.slots() - self.params.ringq().degree()];
        }

        for i in 0..self.params.ringq().degree() {
            let c = fp_ecd[i].round() as i64;
            if c >= 0 {
                for j in 0..self.params.ringq().level() {
                    p_out.coeffs[j][i] = c as u64;
                }
            } else {
                for j in 0..self.params.ringq().level() {
                    p_out.coeffs[j][i] = (c + self.params.ringq().modulus_idx(j) as i64) as u64;
                }
            }
        }

        self.params.ringq().ntt_inplace(p_out);
    }

    /// Encodes a bigint vector to multiple polynomials with randomization and returns it.
    /// The number of polynomials is v.len() / slots.
    pub fn encode_randomize_chunk(&mut self, v: &[Integer], std_dev: f64) -> Vec<Poly> {
        let mut p_out = vec![self.params.ringq().new_poly(); v.len() / self.params.slots()];
        self.encode_randomize_chunk_assign(v, std_dev, &mut p_out);
        return p_out;
    }

    /// Encodes a bigint vector to multiple polynomials with randomization and assigns to p_out.
    /// The number of polynomials is v.len() / slots.
    pub fn encode_randomize_chunk_assign(
        &mut self,
        v: &[Integer],
        std_dev: f64,
        p_out: &mut [Poly],
    ) {
        for i in 0..v.len() / self.params.slots() {
            self.encode_randomize_assign(
                &v[i * self.params.slots()..(i + 1) * self.params.slots()],
                std_dev,
                &mut p_out[i],
            );
        }
    }

    /// Decodes a polynomial to bigint vector and returns it.
    pub fn decode(&self, p: &Poly) -> Vec<Integer> {
        let mut v = vec![Integer::ZERO; self.params.slots()];
        self.decode_assign(p, &mut v);
        return v;
    }

    /// Decodes a polynomial to bigint vector and assigns to v_out.
    pub fn decode_assign(&self, p: &Poly, v_out: &mut [Integer]) {
        let p_big = self.params.ringq().to_bigint_balanced(p);

        for i in 0..self.params.slots() {
            v_out[i].assign(0);
            for j in (0..self.params.digits()).rev() {
                v_out[i] *= &self.base_big;
                v_out[i] += &p_big[j * self.params.slots() + i];
            }
            v_out[i].modulo_mut(&self.params.modulus());
        }
    }

    /// Decodes multiple polynomials to bigint vector and returns it.
    /// The length of output vector is p.len() * slots.
    pub fn decode_chunk(&self, p: &[Poly]) -> Vec<Integer> {
        let mut v = vec![Integer::ZERO; p.len() * self.params.slots()];
        self.decode_chunk_assign(p, &mut v);
        return v;
    }

    /// Decodes multiple polynomials to bigint vector and assigns to v_out.
    /// The length of v_out is p.len() * slots.
    pub fn decode_chunk_assign(&self, p: &[Poly], v_out: &mut [Integer]) {
        for i in 0..p.len() {
            self.decode_assign(
                &p[i],
                &mut v_out[i * self.params.slots()..(i + 1) * self.params.slots()],
            );
        }
    }
}

/// EncoderFixed is a special type of Encoder that computes
/// the randomized encoding with fixed std_dev.
#[derive(Clone)]
pub struct EncoderFixed {
    pub params: Parameters,

    pub gaussian_sampler: TwinCDTSampler,

    pub base_big: Integer,
    pub delta_inv: Vec<f64>,
}

impl EncoderFixed {
    /// Creates a new encoder with fixed std_dev.
    pub fn new(params: &Parameters, std_dev: f64) -> EncoderFixed {
        let base_big = Integer::from(params.modulus_base());

        let prec = params.modulus().significant_bits();
        let base_float: Float = Float::with_val(prec, params.modulus_base());
        let mut modulus_inv: Float = -1 / Float::with_val(prec, params.modulus());
        let mut delta_inv = vec![0.0f64; params.digits()];
        for i in 0..delta_inv.len() {
            delta_inv[i] = modulus_inv.to_f64();
            modulus_inv *= &base_float;
        }

        EncoderFixed {
            params: params.clone(),
            gaussian_sampler: TwinCDTSampler::new(params, std_dev),
            base_big,
            delta_inv,
        }
    }

    /// Encodes a bigint vector to polynomial with randomization and returns it.
    pub fn encode_randomize(&mut self, v: &[Integer]) -> Poly {
        let mut p_out = self.params.ringq().new_poly();
        self.encode_randomize_assign(v, &mut p_out);
        return p_out;
    }

    /// Encodes a bigint vector to polynomial with randomization and assigns to p_out.
    pub fn encode_randomize_assign(&mut self, v: &[Integer], p_out: &mut Poly) {
        let mut fp_ecd = vec![0.0f64; self.params.ringq().degree()];

        let mut coeff = Integer::ZERO;
        let mut rem = Integer::ZERO;
        for i in 0..v.len() {
            coeff.assign(&v[i]);
            coeff.modulo_mut(&self.params.modulus());
            for j in 0..(self.params.digits() - 1) {
                rem.assign(&self.base_big);
                coeff.div_rem_mut(&mut rem);
                let r = rem.to_f64();
                fp_ecd[j * self.params.slots() + i] = r;
            }
            let r = coeff.to_f64();
            fp_ecd[(self.params.digits() - 1) * self.params.slots() + i] = r;
        }

        let mut fp_sample = vec![0.0f64; self.params.ringq().degree()];
        for i in 0..self.params.digits() {
            let d = self.params.ringq().degree() - (i + 1) * self.params.slots();
            for j in 0..self.params.ringq().degree() - d {
                fp_sample[j + d] += self.delta_inv[i] * fp_ecd[j];
            }
            for j in self.params.ringq().degree() - d..self.params.ringq().degree() {
                fp_sample[j + d - self.params.ringq().degree()] -= self.delta_inv[i] * fp_ecd[j];
            }
        }

        for i in 0..self.params.ringq().degree() {
            fp_sample[i] = self.gaussian_sampler.read_coset_f64(fp_sample[i]);
        }

        let base_f64 = self.params.modulus_base() as f64;
        for i in 0..self.params.ringq().degree() - self.params.slots() {
            fp_ecd[i + self.params.slots()] =
                fp_sample[i] - base_f64 * fp_sample[i + self.params.slots()];
        }
        for i in self.params.ringq().degree() - self.params.slots()..self.params.ringq().degree() {
            fp_ecd[i + self.params.slots() - self.params.ringq().degree()] = -fp_sample[i]
                - base_f64 * fp_sample[i + self.params.slots() - self.params.ringq().degree()];
        }

        for i in 0..self.params.ringq().degree() {
            let c = fp_ecd[i].round() as i64;
            if c >= 0 {
                for j in 0..self.params.ringq().level() {
                    p_out.coeffs[j][i] = c as u64;
                }
            } else {
                for j in 0..self.params.ringq().level() {
                    p_out.coeffs[j][i] = (c + self.params.ringq().modulus_idx(j) as i64) as u64;
                }
            }
        }

        self.params.ringq().ntt_inplace(p_out);
    }

    /// Encodes a bigint vector to multiple polynomials with randomization and returns it.
    /// The number of polynomials is v.len() / slots.
    pub fn encode_randomize_chunk(&mut self, v: &[Integer]) -> Vec<Poly> {
        let mut p_out = vec![self.params.ringq().new_poly(); v.len() / self.params.slots()];
        self.encode_randomize_chunk_assign(v, &mut p_out);
        return p_out;
    }

    /// Encodes a bigint vector to multiple polynomials with randomization and assigns to p_out.
    /// The number of polynomials is v.len() / slots.
    pub fn encode_randomize_chunk_assign(&mut self, v: &[Integer], p_out: &mut [Poly]) {
        for i in 0..v.len() / self.params.slots() {
            self.encode_randomize_assign(
                &v[i * self.params.slots()..(i + 1) * self.params.slots()],
                &mut p_out[i],
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::ParametersLiteral;

    use super::*;

    #[test]
    fn test_encode() {
        let params = Parameters::new(ParametersLiteral::logn_19_logp_256());
        let mut oracle = UniformSampler::new(&params);
        let encoder = Encoder::new(&params);

        let mut v = vec![Integer::ZERO; params.slots()];
        for i in 0..params.slots() {
            oracle.read_big_mod_assign(&mut v[i]);
        }

        let p_out = encoder.encode(&v);
        let v_out = encoder.decode(&p_out);

        for i in 0..params.slots() {
            assert_eq!(v[i], v_out[i]);
        }
    }

    #[test]
    fn test_encode_randomize() {
        let params = Parameters::new(ParametersLiteral::logn_19_logp_256());
        let mut oracle = UniformSampler::new(&params);
        let mut encoder = Encoder::new(&params);

        let mut v = vec![Integer::ZERO; params.slots()];
        for i in 0..params.slots() {
            oracle.read_big_mod_assign(&mut v[i]);
        }

        let p_out = encoder.encode_randomize(&v, params.blind_std_dev());
        let v_out = encoder.decode(&p_out);

        for i in 0..params.slots() {
            assert_eq!(v[i], v_out[i]);
        }
    }

    #[test]
    fn test_encode_fixed_randomize() {
        let params = Parameters::new(ParametersLiteral::logn_19_logp_256());
        let mut oracle = UniformSampler::new(&params);
        let encoder = Encoder::new(&params);
        let mut encoder_fixed = EncoderFixed::new(&params, params.commit_std_dev());

        let mut v = vec![Integer::ZERO; params.slots()];
        for i in 0..params.slots() {
            oracle.read_big_mod_assign(&mut v[i]);
        }

        let p_out = encoder_fixed.encode_randomize(&v);
        let v_out = encoder.decode(&p_out);

        for i in 0..params.slots() {
            assert_eq!(v[i], v_out[i]);
        }
    }
}
