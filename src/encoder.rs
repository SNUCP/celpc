use super::csprng::*;
use super::ring::*;
use super::*;
use ethnum::U256;
use rug::integer::Order;
use rug::{ops::*, Assign, Integer};

pub struct Encoder<'a> {
    pub params: &'a Parameters,
    pub sampler: KarneySampler,
}

impl<'a> Encoder<'a> {
    /// Creates a new encoder.
    pub fn new(params: &'a Parameters) -> Encoder<'a> {
        Encoder {
            params: params,
            sampler: KarneySampler::new(),
        }
    }

    /// Encodes a vector of U256 into a polynomial.
    /// We encode s elements of Z_p into a single polynomial in R_q.
    pub fn encode(&self, v: &[U256]) -> Poly {
        let mut pout = self.params.ringq.new_poly();
        self.encode_assign(v, &mut pout);
        return pout;
    }

    /// Encodes a vector of U256 into a polynomial.
    /// We encode s elements of Z_p into a single polynomial in R_q.
    pub fn encode_assign(&self, v: &[U256], pout: &mut Poly) {
        let params = &self.params;
        pout.clear();

        let b256 = U256::from(params.b);
        for (i, a) in v.iter().enumerate() {
            let mut amod = a % self.params.p;
            for j in 0..self.params.kap - 1 {
                let ac = (amod % b256).as_u64();
                for k in 0..params.q.len() {
                    pout.coeffs[k][i + j * params.s] = ac;
                }
                amod /= b256;
            }
            let ac = amod.as_u64();
            for k in 0..params.q.len() {
                pout.coeffs[k][i + params.s * (params.kap - 1)] = ac;
            }
        }

        pout.is_ntt = false;
        params.ringq.ntt(pout);
    }

    /// Encodes a chunk of vectors of U256 into a chunk of polynomials.
    pub fn encode_chunk_assign(&mut self, v: &[U256], pout: &mut [Poly]) {
        if v.len() != pout.len() * self.params.s {
            panic!("invalid length");
        }
        for (i, p) in pout.iter_mut().enumerate() {
            self.encode_assign(&v[i * self.params.s..(i + 1) * self.params.s], p);
        }
    }

    /// Computes pout += p * cx^d.
    /// d must be smaller than p.len().
    #[inline]
    fn monomial_mul_and_add_assign(&self, p: &[f64], c: f64, d: usize, pout: &mut [f64]) {
        let n = p.len();
        for i in 0..n - d {
            pout[i + d] += c * p[i];
        }
        for i in n - d..n {
            pout[i + d - n] -= c * p[i];
        }
    }

    /// Encodes a vector of U256 into a polynomial, with gaussian noise.
    /// We encode s elements of Z_p into a single polynomial in R_q.
    pub fn encode_randomized(&mut self, v: &[U256], sigma: f64) -> Poly {
        let mut pout = self.params.ringq.new_poly();
        self.encode_randomized_assign(v, sigma, &mut pout);
        return pout;
    }

    /// Encodes a vector of U256 into a polynomial, with gaussian noise.
    /// We encode s elements of Z_p into a single polynomial in R_q.
    pub fn encode_randomized_assign(&mut self, v: &[U256], sigma: f64, pout: &mut Poly) {
        let params = self.params;
        pout.clear();

        let mut buff0 = vec![0.0; params.d];
        let mut buff1 = vec![0.0; params.d];

        // Encode v to float
        let bf64 = self.params.b as f64;
        let b256 = U256::from(params.b);
        for (i, a) in v.iter().enumerate() {
            let mut amod = a % self.params.p;
            for j in 0..self.params.kap - 1 {
                buff0[i + j * params.s] = (amod % b256).as_f64();
                amod /= b256;
            }
            buff0[i + params.s * (params.kap - 1)] = amod.as_f64();
        }

        // Multiply P^-1 = -1/(b^d/s + 1) (X^(d-s) + b*X^(d-2s) + b^2 X^(d-3s) + ... + b^(d/s-1))
        let mut pinv = -1.0 / (params.p.as_f64());
        for i in 1..=params.kap {
            self.monomial_mul_and_add_assign(&buff0, pinv, params.d - i * params.s, &mut buff1);
            pinv *= bf64;
        }

        // Sample a* from coset P^-1 * a.
        for i in 0..params.d {
            buff1[i] = self.sampler.sample_coset(buff1[i], sigma);
        }

        // Compute (X^s - b) * a*.
        for i in 0..params.d - params.s {
            buff0[i + params.s] = buff1[i] - bf64 * buff1[i + params.s];
        }
        for i in params.d - params.s..params.d {
            buff0[i + params.s - params.d] = -buff1[i] - bf64 * buff1[i + params.s - params.d];
        }

        // Finally, put result into pOut.
        for i in 0..buff0.len() {
            let c = buff0[i].round() as i64;
            if c < 0 {
                for j in 0..params.q.len() {
                    pout.coeffs[j][i] = (c + (params.q[j] as i64)) as u64;
                }
            } else {
                for j in 0..params.q.len() {
                    pout.coeffs[j][i] = c as u64;
                }
            }
        }

        pout.is_ntt = false;
        params.ringq.ntt(pout);
    }

    /// Encodes a chunk of vectors of U256 into a chunk of polynomials, with gaussian noise.
    pub fn encode_randomized_chunk_assign(&mut self, v: &[U256], sigma: f64, pout: &mut [Poly]) {
        if v.len() != pout.len() * self.params.s {
            panic!("invalid length");
        }
        for (i, p) in pout.iter_mut().enumerate() {
            self.encode_randomized_assign(&v[i * self.params.s..(i + 1) * self.params.s], sigma, p);
        }
    }

    /// Decodes a polynomial into a vector of U256.
    /// Output is always length s.
    pub fn decode(&self, p: &Poly) -> Vec<U256> {
        let mut vout = vec![U256::ZERO; self.params.s];
        self.decode_assign(p, &mut vout);
        return vout;
    }

    /// Decodes a polynomial into a vector of U256.
    /// vout must be of length s.
    pub fn decode_assign(&self, p: &Poly, vout: &mut [U256]) {
        if p.is_ntt {
            panic!("p must be in normal domain");
        }

        let mut tmp = Integer::from(0);
        for i in 0..self.params.s {
            vout[i] = U256::ZERO;
            tmp.assign(Integer::ZERO);
            for j in (0..self.params.kap).rev() {
                tmp *= self.params.b;
                tmp += self
                    .params
                    .ringq
                    .get_coeff_balanced(p, i + j * self.params.s)
                    .as_i128();
            }
            tmp.rem_euc_assign(&self.params.pbig);
            for (k, &v) in tmp.to_digits::<u128>(Order::LsfLe).iter().enumerate() {
                vout[i].0[k] = v;
            }
        }
    }

    /// Decodes a chunk of polynomials.
    /// vout must be of length s * p.len().
    pub fn decode_chunk_assign(&self, p: &[Poly], vout: &mut [U256]) {
        if vout.len() != p.len() * self.params.s {
            panic!("invalid length");
        }
        for (i, p) in p.iter().enumerate() {
            self.decode_assign(p, &mut vout[i * self.params.s..(i + 1) * self.params.s]);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::csprng::*;
    use crate::*;
    use ethnum::U256;

    #[test]
    pub fn test_encoder() {
        let params = Parameters::N_19();
        let mut ecd = Encoder::new(&params);

        let mut us = UniformSampler::new();

        let mut msg = vec![U256::ZERO; params.s];
        for j in 0..params.s {
            msg[j] = us.sample_range_u256(params.p);
        }
        let mut m = ecd.encode(&msg);
        params.ringq.intt(&mut m);
        let mout = ecd.decode(&m);
        assert_eq!(msg, mout);

        let mut mr = ecd.encode_randomized(&msg, params.s3);
        params.ringq.intt(&mut mr);
        let mout = ecd.decode(&mr);
        assert_eq!(msg, mout);
    }
}
