use std::f64::consts::PI;

use super::UniformSampler;
use crate::ring::*;
use rug::{
    ops::{Pow, PowAssign},
    Float,
};

/// CDTSampler is a Gaussian sampler for small, fixed parameters.
/// Additionally, it is used as a base sampler for ConvSampler.
///
/// The target precision of the output is 64 bits.
pub struct CDTSampler {
    pub base_sampler: UniformSampler,

    pub c: f64,
    pub s: f64,

    // Probabilities are stored with scaling factor 2^64.
    pub table: Vec<u64>,

    pub tailcut_low: i64,
    pub tailcut_high: i64,
}

impl CDTSampler {
    pub fn new(center: f64, sigma: f64) -> CDTSampler {
        let base_sampler = UniformSampler::new();

        let tailcut_low = (center - 5.0 * sigma).round() as i64;
        let tailcut_high = (center + 5.0 * sigma).round() as i64;
        let table_size = (tailcut_high - tailcut_low + 1) as usize;

        let c_big = Float::with_val(64, center);
        let s_big = Float::with_val(64, sigma);
        let s_big_sq_inv = 1 / s_big.pow(2);

        let mut gaussian_sum = Float::new(64);
        let mut gaussians = Vec::new();
        for x in tailcut_low..tailcut_high {
            let mut rho = Float::with_val(64, x);
            rho -= &c_big;
            rho.pow_assign(2);
            rho *= &s_big_sq_inv;
            rho *= -PI;
            rho.exp_mut();

            gaussians.push(rho);
            gaussian_sum += &gaussians[(x - tailcut_low) as usize];
        }

        let exp64 = Float::with_val(64, 1u128 << 64);
        for i in 0..(table_size - 1) {
            gaussians[i] /= &gaussian_sum;
            gaussians[i] *= &exp64;
        }

        let mut table = vec![0u64; table_size];
        table[table_size - 1] = u64::MAX;
        gaussian_sum = gaussians[0].clone();
        for i in 1..(table_size - 1) {
            table[i] = gaussian_sum
                .to_integer()
                .unwrap()
                .to_u64()
                .unwrap_or(u64::MAX);
            gaussian_sum += &gaussians[i];
        }

        CDTSampler {
            base_sampler: base_sampler,
            c: center,
            s: sigma,

            table: table,
            tailcut_low,
            tailcut_high,
        }
    }

    /// Returns a sample from the Gaussian distribution.
    pub fn sample(&mut self) -> i64 {
        let r = self.base_sampler.sample_u64();

        let x = match self.table.binary_search(&r) {
            Ok(x) => x,
            Err(x) => x - 1,
        };
        return x as i64 + self.tailcut_low;
    }

    /// Samples a polynomial with coefficients from the Gaussian distribution.
    /// Output polynomial is in NTT domain.
    pub fn sample_poly_assign(&mut self, ring: &Ring, pout: &mut Poly) {
        pout.is_ntt = false;

        for i in 0..ring.degree {
            let cc = self.sample();
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

/// CDTSamplerVarCenter is a variant of CDTSampler with variable center and fixed sigma.
/// This is possible by using SampleC algorithm from [MW18].
/// Currently, base = 32 and precision = 6, which implies statistical distance ~ 2^-60.
pub struct CDTSamplerVarCenter {
    pub base_samplers: Vec<CDTSampler>,

    pub b_log: u64,
    pub k: u64,
    pub prec_log: u64,
}

impl CDTSamplerVarCenter {
    pub fn new(sigma: f64) -> CDTSamplerVarCenter {
        let b_log = 5;
        let k = 6;

        let mut base_samplers = Vec::new();
        for i in 0..(1 << b_log) {
            base_samplers.push(CDTSampler::new((i as f64) / f64::exp2(b_log as f64), sigma));
        }

        CDTSamplerVarCenter {
            base_samplers: base_samplers,

            b_log: b_log,
            k: k,
            prec_log: b_log * k,
        }
    }

    /// Returns a sample from the Gaussian distribution.
    pub fn sample(&mut self, c: f64) -> i64 {
        let ci = c.floor() as i64;
        let c = c - c.floor(); // c is now always in [0, 1).

        let f64_prec_log = 52;
        let tail_prec_log = f64_prec_log - self.prec_log;

        // Scale 53 bits...
        let cx = (c * f64::exp2(f64_prec_log as f64)) as u64;

        // Extract first 30 bits of precision
        let cmsb = cx >> tail_prec_log;
        // And the last 22 bits.
        let clsb = cx & ((1 << tail_prec_log) - 1);

        let mut x = self.sampleC(cmsb);

        // Use Bernoulli to sample the rounding.
        let r = self.base_samplers[0].base_sampler.sample_u64() & ((1 << tail_prec_log) - 1);
        if clsb >= r {
            x += 1;
        }

        return x + ci;
    }

    /// Samples from Gaussian distribution with center c / 2^30.
    fn sampleC(&mut self, c: u64) -> i64 {
        let mut x = c as i64;
        let mask = (1 << self.b_log) - 1;
        for _ in 0..self.k {
            let mut g = self.base_samplers[(x & mask) as usize].sample();
            if (x & mask) > 0 && x < 0 {
                g -= 1;
            }
            x >>= self.b_log;
            x += g;
        }
        return x;
    }

    /// Returns a sample from the Gaussian distribution over a coset c + Z.
    pub fn sample_coset(&mut self, center: f64) -> f64 {
        return center + (self.sample(-center) as f64);
    }
}
