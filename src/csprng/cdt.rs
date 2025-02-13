use core::f64;
use std::f64::consts::PI;

use super::UniformSampler;
use crate::ring::*;
use rug::{ops::PowAssign, Assign, Float};

const TAIL_CUT: f64 = 6.0;

pub fn generate_cdt(center: f64, sigma: f64) -> Vec<u64> {
    let c_big = Float::with_val(64, center);
    let s_big = Float::with_val(64, sigma);
    let s_big_sq_inv = 1 / Float::with_val(64, sigma * sigma);

    let tailcut_high = (TAIL_CUT * sigma).ceil() as i64;
    let tailcut_low = -tailcut_high;
    let table_size = (tailcut_high - tailcut_low + 1) as usize;

    let mut table = vec![0u64; table_size];
    let mut cdf = Float::with_val(64, 0);
    let mut cdf_buf = Float::with_val(64, 0);
    let exp_64 = Float::with_val(64, 1u128 << 64);
    for x in tailcut_low..=tailcut_high {
        let i = (x - tailcut_low) as usize;

        let mut rho = Float::with_val(64, x);
        rho -= &c_big;
        rho.pow_assign(2);
        rho *= &s_big_sq_inv;
        rho *= -PI;
        rho.exp_mut();
        rho /= &s_big;

        cdf += &rho;

        cdf_buf.assign(&cdf);
        cdf_buf *= &exp_64;
        table[i] = cdf_buf.to_integer().unwrap().to_u64().unwrap_or(u64::MAX);
    }

    return table;
}

/// CDTSampler is a Gaussian sampler for small, fixed parameters.
/// Additionally, it is used as a base sampler for ConvSampler.
///
/// The target precision of the output is 64 bits.
pub struct CDTSampler {
    pub base_sampler: UniformSampler,

    pub center: f64,
    pub sigma: f64,

    table: Vec<u64>,
    center_floor: i64,
    tailcut_low: i64,
}

impl CDTSampler {
    pub fn new(center: f64, sigma: f64) -> CDTSampler {
        let base_sampler = UniformSampler::new();

        let tailcut_low = -(TAIL_CUT * sigma).ceil() as i64;
        let center_floor = center.floor();
        let table = generate_cdt(center - center_floor, sigma);

        CDTSampler {
            base_sampler: base_sampler,
            center,
            sigma,

            table,
            center_floor: center_floor as i64,
            tailcut_low,
        }
    }

    /// Returns a sample from the Gaussian distribution.
    pub fn sample(&mut self) -> i64 {
        let r = self.base_sampler.sample_u64();
        let x = self.table.binary_search(&r).unwrap_or_else(|x| x - 1);
        return x as i64 + self.center_floor + self.tailcut_low;
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
/// Currently, base = 256 and precision = 4, which implies statistical distance ~ 2^-64.
pub struct CDTSamplerVarCenter {
    pub base_samplers: Vec<CDTSampler>,

    pub b_log: u64,
    pub k: u64,
    pub prec_log: u64,
}

impl CDTSamplerVarCenter {
    pub fn new(sigma: f64) -> CDTSamplerVarCenter {
        let b_log = 8;
        let k = 4;

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

        // Extract first bits of precision
        let cmsb = cx >> tail_prec_log;
        // And the last bits.
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

/// TwinCDTSampler is a variant of CDTSampler with variable center and fixed sigma, based on [AR18].
/// Currently, base = 128.
pub struct TwinCDTSampler {
    pub base_sampler: UniformSampler,

    pub sigma: f64,

    base: usize,
    tables: Vec<Vec<u64>>,
    tailcut_low: i64,
}

impl TwinCDTSampler {
    pub fn new(sigma: f64) -> TwinCDTSampler {
        let base = 256;

        let tailcut_low = -(TAIL_CUT * sigma).ceil() as i64;
        let mut tables = Vec::with_capacity(base);
        for i in 0..base {
            tables.push(generate_cdt((i as f64) / (base as f64), sigma));
        }

        return TwinCDTSampler {
            base_sampler: UniformSampler::new(),

            sigma,

            base,
            tables,
            tailcut_low,
        };
    }

    pub fn sample(&mut self, center: f64) -> i64 {
        let c_floor = center.floor();

        let c0 = ((self.base as f64) * (center - c_floor)).floor() as usize % self.base;
        let c1 = ((self.base as f64) * (center - c_floor)).ceil() as usize % self.base;

        let u = self.base_sampler.sample_u64();

        let v0 = self.tables[c0].binary_search(&u).unwrap_or_else(|x| x - 1);
        let v1 = self.tables[c1].binary_search(&u).unwrap_or_else(|x| x - 1);

        if v0 == v1 {
            return v0 as i64 + self.tailcut_low + c_floor as i64;
        }

        let c_frac = center - c_floor;
        let mut cdf = 0.0;
        for x in self.tailcut_low..=(v0 as i64) {
            let xf = x as f64;
            cdf += f64::exp(-PI * (xf - c_frac) * (xf - c_frac) / (self.sigma * self.sigma))
                / self.sigma;
        }

        let p = (u as f64) / f64::exp2(64.0);
        if p < cdf {
            return v0 as i64 + self.tailcut_low + c_floor as i64;
        }
        return v1 as i64 + self.tailcut_low + c_floor as i64;
    }

    /// Returns a sample from the Gaussian distribution over a coset c + Z.
    pub fn sample_coset(&mut self, center: f64) -> f64 {
        return center + (self.sample(-center) as f64);
    }
}
