use std::f64::consts::SQRT_2;

use rug::{ops::Pow, Float};

use super::CDTSamplerVarCenter;
use crate::ring::*;

/// ConvolveSampler is a Gaussian sampler for large, variable parameters.
/// It is based on [MW18].
///
/// The target precision of the output is 52 bits.
pub struct ConvolveSampler {
    pub base_sampler: CDTSamplerVarCenter,

    pub sigs: Vec<f64>,
    pub weights: Vec<i64>,

    pub K: f64,
}

impl ConvolveSampler {
    pub fn new(sigma: f64) -> ConvolveSampler {
        let eta = 6.0;
        let mut sigs = Vec::new();
        let mut weights = Vec::new();

        sigs.push(34.0);
        weights.push(0);

        let factor = 3;

        for i in 1..factor {
            let z = (sigs[i - 1] / (SQRT_2 * eta)).floor() as i64;
            let s =
                ((((z * z) + i64::max((z - 1) * (z - 1), 1)) as f64) * sigs[i - 1] * sigs[i - 1])
                    .sqrt();

            sigs.push(s);
            weights.push(z);
        }

        let base_sampler = CDTSamplerVarCenter::new(sigs[0]);

        let b_big = Float::with_val(64, 1 << base_sampler.b_log);
        let mut sbar_big = Float::with_val(64, 0);
        for i in 0..(base_sampler.k as i64) {
            sbar_big += b_big.clone().pow(-2 * i);
        }
        sbar_big.sqrt_mut();
        sbar_big *= sigs[0];

        let s_big = Float::with_val(64, sigma);
        let mut K_big: Float = s_big.pow(2) - sbar_big.pow(2);
        K_big.sqrt_mut();
        K_big /= sigs[factor - 1];

        ConvolveSampler {
            base_sampler: base_sampler,

            sigs: sigs,
            weights: weights,

            K: K_big.to_f64(),
        }
    }

    /// Algorithm SampleI from [MW18].
    pub fn sampleI(&mut self, i: usize) -> i64 {
        if i == 0 {
            return self.base_sampler.base_samplers[0].sample();
        }

        let x1 = self.sampleI(i - 1);
        let x2 = self.sampleI(i - 1);

        let z1 = self.weights[i];
        let z2 = i64::max(1, z1 - 1);

        return z1 * x1 + z2 * x2;
    }

    /// Returns a sample from the Gaussian distribution.
    pub fn sample(&mut self, c: f64) -> i64 {
        let x = self.sampleI(self.sigs.len() - 1);

        let cc = c + self.K * (x as f64);

        return self.base_sampler.sample(cc);
    }

    /// Returns a sample from the Gaussian distribution over a coset c + Z.
    pub fn sample_coset(&mut self, c: f64) -> f64 {
        return c + (self.sample(-c) as f64);
    }

    /// Samples a polynomial with coefficients from the Gaussian distribution.
    /// Output polynomial is in NTT domain.
    pub fn sample_poly_assign(&mut self, ring: &Ring, pout: &mut Poly) {
        pout.is_ntt = false;

        for i in 0..ring.degree {
            let cc = self.sample(0.0);
            if cc < 0 {
                for j in 0..ring.moduli.len() {
                    pout.coeffs[j][i] = (cc + (ring.moduli[j] as i64)) as u64;
                }
            } else {
                for j in 0..ring.moduli.len() {
                    pout.coeffs[j][i] = cc as u64;
                }
            }
        }
        ring.ntt(pout);
    }
}
