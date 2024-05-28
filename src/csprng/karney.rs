use std::f64::consts::*;

use crate::ring::*;

use super::UniformSampler;

// Helpers from https://eprint.iacr.org/2019/068

/// returns a+b <= c.
#[inline]
fn is_a_plus_b_leq_c(a: f64, b: f64, c: f64) -> bool {
    if a > c || b > c {
        return false;
    }
    if a > c / 2.0 {
        return b <= c - a;
    }
    if b > c / 2.0 {
        return a <= c - b;
    }
    return true;
}

/// computes a+b >= c.
#[inline]
fn is_a_plus_b_geq_c(a: f64, b: f64, c: f64) -> bool {
    if a < c || b < c {
        return false;
    }
    if a < c / 2.0 {
        return b >= c - a;
    }
    if b < c / 2.0 {
        return a >= c - b;
    }
    return true;
}

/// returns floor(ta) and (ta) % 1.
#[inline]
fn mul_floor_fract(t: i64, a: f64) -> (i64, f64) {
    let s = t.signum();
    let mut tu = t.abs() as u64;

    let af = a.fract();
    let mut ta = 0.0;
    let mut taf = 0.0;
    let mut i = 0;
    while tu > 0 {
        if tu & 1 == 1 {
            ta += float_extras::f64::ldexp(a, i);
            taf += float_extras::f64::ldexp(af, i);
        }
        tu >>= 1;
        i += 1;
    }
    if s < 0 {
        ta = -ta;
        taf = -taf;
    }

    let x0 = t * (a.floor() as i64) + taf.floor() as i64;
    let x1 = ta.fract();
    return (x0, x1);
}

#[inline]
fn compute_i(sig: f64, c: f64, t: i64, s: i64) -> i64 {
    let (mut i, b) = mul_floor_fract(t, sig);
    if s < 0 && b > c {
        i += 1;
    }
    if s > 0 && (b > 0.0 || c > 0.0) {
        if is_a_plus_b_leq_c(b, c, 1.0) {
            i += 1;
        } else {
            i += 2;
        }
    }
    return i;
}

#[inline]
fn is_x_bar_zero(sig: f64, c: f64, t: i64, s: i64) -> bool {
    let (_, b) = mul_floor_fract(t, sig);
    if s > 0 {
        if (b == 0.0 && c == 0.0) || (b + c == 1.0) {
            return true;
        }
    }
    if s < 0 && b == c {
        return true;
    }
    return false;
}

#[inline]
fn check_reject_a(sig: f64, c: f64, t: i64, s: i64, j: i64) -> bool {
    if (j < sig.floor() as i64) || is_x_bar_zero(sig, c, t, s) {
        return false;
    }

    let (_, b) = mul_floor_fract(t, sig);
    let a = sig.fract();
    let z = a + b;
    if s > 0 {
        if is_a_plus_b_geq_c(b, c, 1.0) {
            return is_a_plus_b_leq_c(z, c, 2.0);
        } else {
            return is_a_plus_b_leq_c(z, c, 1.0);
        }
    } else {
        if b > c {
            if z <= 1.0 {
                return true;
            } else {
                return c >= z.fract();
            }
        } else {
            return c >= z;
        }
    }
}

#[inline]
fn check_reject_b(c: f64, t: i64, s: i64, j: i64) -> bool {
    return c == 0.0 && t == 0 && j == 0 && s < 0;
}

#[inline]
fn compute_x(sig: f64, c: f64, t: i64, s: i64, j: i64) -> f64 {
    if is_x_bar_zero(sig, c, t, s) {
        return (j as f64) / sig;
    }

    let (_, b) = mul_floor_fract(t, sig);
    let mut xbar = 1.0 - (b + (s as f64) * c);
    if s < 0 && b < c {
        xbar -= 1.0;
    }
    if s > 0 && is_a_plus_b_geq_c(b, c, 1.0) {
        xbar += 1.0;
    }
    return (xbar + (j as f64)) / sig;
}

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
            let t = self.normal_sample();
            let s = 2 * (self.base_sampler.sample_i64() & 1) - 1;
            let j = self.base_sampler.sample_range(sigma.ceil() as u64) as i64;

            if check_reject_a(sigma, c, t, s, j) {
                continue;
            }
            if check_reject_b(c, t, s, j) {
                continue;
            }

            let i = compute_i(sigma, c, t, s);
            let x = compute_x(sigma, c, t, s, j);

            let p = f64::exp(-0.5 * x * (2.0 * (t as f64) + x));
            if self.base_sampler.sample_f64() <= p {
                return s * ((i as i64) + j);
            }
        }
    }

    /// Samples a random i64 with integer mean c and standard deviation sigma.
    /// This sets sigma / sqrt(2 * PI) up to the nearest integer,
    /// so the resulting sample may have slightly larger standard deviation.
    /// However, this is faster than normal sample.
    pub fn sample_i64_exact(&mut self, c: i64, sigma: f64) -> i64 {
        let sigma = (sigma as f64 / f64::sqrt(2.0 * PI)).ceil() as i64;
        loop {
            let t = self.normal_sample();
            let s = 2 * (self.base_sampler.sample_i64() & 1) - 1;

            let i = t * sigma + s * c; // always integer, so x = j / sigma
            let j = self.base_sampler.sample_range(sigma as u64) as i64;
            let x = (j, sigma);
            if x.0 >= x.1 {
                continue;
            }

            if t == 0 && x.0 == 0 && s < 0 {
                continue;
            }

            let x = (x.0 as f64) / (x.1 as f64);
            let p = f64::exp(-0.5 * x * (2.0 * (t as f64) + x));
            if self.base_sampler.sample_f64() <= p {
                return s * (i + j);
            }
        }
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

        for i in 0..ring.degree {
            let cc = self.sample_i64(c, sigma);
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

    /// Samples a random polynomial. See sample_i64_exact.
    /// Output polynomial is in NTT domain.
    pub fn sample_poly_exact_assign(&mut self, ring: &Ring, c: i64, sigma: f64, pout: &mut Poly) {
        pout.is_ntt = false;

        for i in 0..ring.degree {
            let cc = self.sample_i64_exact(c, sigma);
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
