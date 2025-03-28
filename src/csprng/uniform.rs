use std::{
    cmp::Ordering,
    io::{Read, Write},
};

use rand::prelude::*;
use rug::{integer::Order, Integer};
use sha3::{
    digest::{ExtendableOutput, ExtendableOutputReset, Update},
    Shake128, Shake128Reader,
};

use crate::parameters::Parameters;
use crate::ring::*;

/// BUFSIZE is the default buffer size of UniformSampler.
const BUFSIZE: usize = 8192;

/// UniformSampler samples values from a uniform distribution.
#[derive(Clone)]
pub struct UniformSampler {
    pub params: Parameters,

    pub hasher: Shake128,
    pub xof: Shake128Reader,

    pub buf: [u8; BUFSIZE],
    pub ptr: usize,

    mod_max: Integer,
    mod_buf: Vec<u8>,
    msb_mask: u8,
}

impl UniformSampler {
    /// Creates a new uniform sampler.
    pub fn new(params: &Parameters) -> UniformSampler {
        let mut seed = [0u8; 32];

        let mut rng = StdRng::from_entropy();
        rng.fill_bytes(&mut seed);

        UniformSampler::new_with_seed(params, &seed)
    }

    /// Creates a new uniform sampler with seed.
    pub fn new_with_seed(params: &Parameters, seed: &[u8]) -> UniformSampler {
        let mut hasher = Shake128::default();
        hasher.update(seed);
        let xof = hasher.clone().finalize_xof();

        let mod_max = params.modulus().clone() - Integer::ONE;
        let k = (mod_max.significant_bits() + 7) / 8;
        let b = mod_max.significant_bits() % 8;

        let mod_buf = vec![0u8; k as usize];
        let msb_mask = (1 << b) - 1;

        UniformSampler {
            params: params.clone(),

            hasher: hasher,
            xof: xof,

            buf: [0u8; BUFSIZE],
            ptr: BUFSIZE,

            mod_max: mod_max,
            mod_buf: mod_buf,
            msb_mask: msb_mask,
        }
    }

    /// Creates a new UniformSampler as random oracle.
    pub fn new_oracle(params: &Parameters) -> UniformSampler {
        UniformSampler::new_with_seed(params, &[])
    }

    /// reset resets the internal state of the sampler.
    pub fn reset(&mut self) {
        self.hasher.finalize_xof_reset();
        self.xof = self.hasher.clone().finalize_xof();
        self.ptr = BUFSIZE;
    }

    /// finalize finalizes the internal state to XOF.
    pub fn finalize(&mut self) {
        self.xof = self.hasher.clone().finalize_xof();
        self.ptr = BUFSIZE;
    }

    /// write_poly writes a polynomial to the sampler.
    pub fn write_poly(&mut self, p: &Poly) {
        for i in 0..self.params.ringq().level() {
            for j in 0..self.params.ringq().degree() {
                self.hasher.write(&p.coeffs[i][j].to_le_bytes()).unwrap();
            }
        }
    }

    /// read_u64 generates a random u64.
    pub fn read_u64(&mut self) -> u64 {
        if self.ptr == BUFSIZE {
            self.xof.read_exact(&mut self.buf).unwrap();
            self.ptr = 0;
        }

        let mut res = 0;
        res |= (self.buf[self.ptr + 0] as u64) << (0 * 8);
        res |= (self.buf[self.ptr + 1] as u64) << (1 * 8);
        res |= (self.buf[self.ptr + 2] as u64) << (2 * 8);
        res |= (self.buf[self.ptr + 3] as u64) << (3 * 8);
        res |= (self.buf[self.ptr + 4] as u64) << (4 * 8);
        res |= (self.buf[self.ptr + 5] as u64) << (5 * 8);
        res |= (self.buf[self.ptr + 6] as u64) << (6 * 8);
        res |= (self.buf[self.ptr + 7] as u64) << (7 * 8);

        self.ptr += 8;

        return res;
    }

    /// read_u64_range generates a random u64 in the range [0, bound).
    pub fn read_u64_range(&mut self, bound: u64) -> u64 {
        let max = u64::MAX - u64::MAX % bound;
        loop {
            let r = self.read_u64();
            if r < max {
                return r % bound;
            }
        }
    }

    /// read_big_mod generates a random big integer in the range [0, modulus).
    pub fn read_big_mod_assign(&mut self, x_out: &mut Integer) {
        loop {
            self.xof.read_exact(&mut self.mod_buf).unwrap();

            self.mod_buf[0] &= self.msb_mask;

            x_out.assign_digits(&self.mod_buf, Order::Msf);
            if (*x_out).cmp(&self.mod_max) == Ordering::Less {
                return;
            }
        }
    }

    /// read_poly_assign generates a random polynomial in NTT form.
    pub fn read_poly_assign(&mut self, p_out: &mut Poly) {
        for i in 0..self.params.ringq().level() {
            let qi = self.params.ringq().modulus_idx(i);
            for j in 0..self.params.ringq().degree() {
                p_out.coeffs[i][j] = self.read_u64_range(qi);
            }
        }
    }

    /// read_monomial_assign generates a random monomial in NTT form.
    pub fn read_monomial_assign(&mut self, p_out: &mut Poly) {
        p_out.clear();

        let d = self.read_u64_range(2 * self.params.ringq().degree() as u64) as usize;

        for i in 0..self.params.ringq().level() {
            if d & 1 == 0 {
                p_out.coeffs[i][d >> 1] = 1;
            } else {
                p_out.coeffs[i][d >> 1] = self.params.ringq().modulus_idx(i) - 1;
            }
        }

        self.params.ringq().ntt_inplace(p_out);
    }
}
