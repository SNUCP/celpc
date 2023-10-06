use std::f64::consts::PI;

use crate::ring::*;

use super::UniformSampler;

pub struct KarneySampler {
    pub base_sampler: UniformSampler,

    // invSqrtE is 1/sqrt(e) = exp(-0.5).
    pub inv_sqrt_e: f64,
}

impl KarneySampler {
    /// Create a new Karney sampler.
    pub fn new() -> KarneySampler {
        KarneySampler {
            base_sampler: UniformSampler::new(),
            inv_sqrt_e: f64::exp(-0.5),
        }
    }

    /// normal_sample samples k with density exp(-k^2/2).
    fn normal_sample(&mut self) -> i64 {
        let mut x = 0;
        let mut will_reject = true;

        while will_reject {
            x = 0;
            will_reject = false;

            while self.base_sampler.sample_f64() < self.inv_sqrt_e {
                x += 1;
            }

            if x < 2 {
                return x;
            }

            for _ in 0..(x * (x - 1)) {
                if self.base_sampler.sample_f64() >= self.inv_sqrt_e {
                    will_reject = true;
                    break;
                }
            }
        }

        return x;
    }

    /// Samples a random i64 with mean c and standard deviation sigma.
    pub fn sample_i64(&mut self, c: f64, sigma: f64) -> i64 {
        let sigma = sigma / f64::sqrt(2.0 * PI);
        loop {
            let k = self.normal_sample();
            let s = 2 * (self.base_sampler.sample_i64() & 1) - 1;

            let tmp = (k as f64) * sigma + (s as f64) * c;
            let i0 = tmp.ceil() as i64;
            let x0 = (tmp.ceil() - tmp) / sigma;
            let j = self.base_sampler.sample_range(sigma.ceil() as u64) as i64;

            let x = x0 + (j as f64) / sigma;
            if x >= 1.0 {
                continue;
            }

            if x == 0.0 {
                if k == 0 && s < 0 {
                    continue;
                } else {
                    return s * (i0 + j);
                }
            }

            let bias = f64::exp(-0.5 * x * (2.0 * (k as f64) + x));
            if self.base_sampler.sample_f64() <= bias {
                return s * (i0 + j);
            }
        }
    }

    /// Samples a random u64 mod q.
    pub fn sample_u64_mod(&mut self, c: f64, sigma: f64, q: u64) -> u64 {
        let x = self.sample_i64(c, sigma);
        x.rem_euclid(q as i64) as u64
    }

    /// Samples a random number from coset c + Z.
    pub fn sample_coset(&mut self, c: f64, sigma: f64) -> f64 {
        return self.sample_i64(-c, sigma) as f64 + c;
    }

    /// Samples a random polynomial.
    /// Output polynomial is in NTT domain.
    pub fn sample_poly(&mut self, ring: &Ring, c: f64, sigma: f64) -> Poly {
        let mut pout = ring.new_poly();
        self.sample_poly_assign(ring, c, sigma, &mut pout);
        return pout;
    }

    /// Samples a random polynomial.
    /// Output polynomial is in NTT domain.
    pub fn sample_poly_assign(&mut self, ring: &Ring, c: f64, sigma: f64, pout: &mut Poly) {
        pout.is_ntt = false;

        for i in 0..pout.coeffs.len() {
            for j in 0..pout.coeffs[i].len() {
                pout.coeffs[i][j] = self.sample_u64_mod(c, sigma, ring.moduli[i]);
            }
        }

        ring.ntt(pout);
    }
}
