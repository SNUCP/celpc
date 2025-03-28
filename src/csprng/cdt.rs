use super::UniformSampler;
use crate::{parameters::Parameters, ring::Poly};

use core::f64;
use core::f64::consts::PI;
const TAIL_CUT: f64 = 9.0;

pub fn generate_cdt(center: f64, sigma: f64) -> Vec<u64> {
    let norm = f64::sqrt(2.0 * PI) * sigma;

    let tailcut_high = (TAIL_CUT * sigma).ceil() as i64;
    let tailcut_low = -tailcut_high;
    let table_size = (tailcut_high - tailcut_low + 1) as usize;

    let mut table = vec![0u64; table_size];
    let mut cdf = 0.0;
    for i in 0..table_size {
        let x = tailcut_low as f64 + i as f64;
        let rho = f64::exp(-(x - center) * (x - center) / (2.0 * sigma * sigma)) / norm;
        cdf += rho;
        if cdf > 1.0 {
            table[i] = u64::MAX;
        } else {
            table[i] = (cdf * f64::exp2(64.0)).round() as u64;
        }
    }

    return table;
}

/// TwinCDTSampler is a variant of CDTSampler with variable center and fixed sigma, based on [AR18].
/// Currently, base = 256.
#[derive(Clone)]
pub struct TwinCDTSampler {
    pub base_sampler: UniformSampler,

    pub std_dev: f64,

    base: usize,
    tables: Vec<Vec<u64>>,
    tailcut_low: i64,
}

impl TwinCDTSampler {
    pub fn new(params: &Parameters, std_dev: f64) -> TwinCDTSampler {
        let base = 256;

        let tailcut_low = -(TAIL_CUT * std_dev).ceil() as i64;
        let mut tables = Vec::with_capacity(base);
        for i in 0..base {
            tables.push(generate_cdt((i as f64) / (base as f64), std_dev));
        }

        return TwinCDTSampler {
            base_sampler: UniformSampler::new(params),
            std_dev,

            base,
            tables,
            tailcut_low,
        };
    }

    /// read_f64 samples random f64 from Discrete Gaussian Distribution.
    pub fn read_i64(&mut self, center: f64) -> i64 {
        let c_floor = center.floor();

        let c0 = ((self.base as f64) * (center - c_floor)).floor() as usize % self.base;
        let c1 = ((self.base as f64) * (center - c_floor)).ceil() as usize % self.base;

        let u = self.base_sampler.read_u64();

        let v0 = self.tables[c0].binary_search(&u).unwrap_or_else(|x| x - 1);
        let v1 = self.tables[c1].binary_search(&u).unwrap_or_else(|x| x - 1);

        if v0 == v1 {
            return v0 as i64 + self.tailcut_low + c_floor as i64;
        }

        let c_frac = center - c_floor;
        let mut cdf = 0.0;
        let norm = f64::sqrt(2.0 * f64::consts::PI) * self.std_dev;
        for x in self.tailcut_low..=(v0 as i64) {
            let xf = x as f64;
            cdf += f64::exp(-(xf - c_frac) * (xf - c_frac) / (2.0 * self.std_dev * self.std_dev))
                / norm;
        }

        let p = (u as f64) / f64::exp2(64.0);
        if p < cdf {
            return v0 as i64 + self.tailcut_low + c_floor as i64;
        }
        return v1 as i64 + self.tailcut_low + c_floor as i64;
    }

    /// read_coset_f64 samples random f64 from Discrete Gaussian Distribution over a coset.
    pub fn read_coset_f64(&mut self, coset: f64) -> f64 {
        return coset + self.read_i64(-coset) as f64;
    }

    /// read_poly_assign samples random f64 from Discrete Gaussian Distribution and assigns to p_out.
    pub fn read_poly_assign(&mut self, center: f64, p_out: &mut Poly) {
        for i in 0..self.base_sampler.params.ringq().degree() {
            let c = self.read_i64(center);
            if c >= 0 {
                for j in 0..self.base_sampler.params.ringq().level() {
                    p_out.coeffs[j][i] = c as u64;
                }
            } else {
                for j in 0..self.base_sampler.params.ringq().level() {
                    p_out.coeffs[j][i] =
                        (c + self.base_sampler.params.ringq().modulus_idx(j) as i64) as u64;
                }
            }
        }

        self.base_sampler.params.ringq().ntt_inplace(p_out);
    }
}
