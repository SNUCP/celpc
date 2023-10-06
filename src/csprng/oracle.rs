use crate::ring::*;
use primitive_types::U256;
use sha3::{
    digest::{ExtendableOutput, Reset, Update, XofReader},
    Shake128, Shake128Reader,
};

pub struct Oracle {
    hasher: Shake128,
    xof: Shake128Reader,
}

impl Oracle {
    /// Create a new uniform sampler.
    pub fn new() -> Oracle {
        Oracle {
            hasher: Shake128::default(),
            xof: Shake128::default().finalize_xof(),
        }
    }

    /// Writes a u256 value to random oracle.
    pub fn write_u256(&mut self, x: U256) {
        let mut buff = [0u8; 32];
        x.to_big_endian(&mut buff);
        self.hasher.update(&buff);
    }

    /// Writes a polynomial to random oracle.
    pub fn write_poly(&mut self, p: &Poly) {
        for i in 0..p.coeffs.len() {
            for j in 0..p.coeffs[i].len() {
                self.hasher.update(&p.coeffs[i][j].to_be_bytes());
            }
        }
    }

    /// Finalizes the random oracle.
    pub fn finalize(&mut self) {
        self.xof = self.hasher.clone().finalize_xof();
        self.hasher.reset();
    }

    /// Reads a u128 value smaller than n.
    /// Return value is u256.
    pub fn read_range(&mut self, n: U256) -> U256 {
        let mut buff = [0u8; 32];
        let bound = U256::MAX - (U256::MAX % n);
        loop {
            self.xof.read(&mut buff);
            let x = U256::from_big_endian(&buff);
            if x < bound {
                return x % n;
            }
        }
    }

    /// Reads a u256 value.
    pub fn read_u256(&mut self) -> U256 {
        let mut buff = [0u8; 32];
        self.xof.read(&mut buff);
        U256::from_big_endian(&buff)
    }

    /// Reads a polynomial.
    /// Output polynomial is in NTT domain.
    pub fn read_poly(&mut self, ring: &Ring) -> Poly {
        let mut pout = ring.new_ntt_poly();
        self.read_poly_assign(ring, &mut pout);
        return pout;
    }

    /// Reads a polynomial.
    pub fn read_poly_assign(&mut self, ring: &Ring, pout: &mut Poly) {
        let mut buff = [0u8; 8];
        for i in 0..pout.coeffs.len() {
            let q = ring.moduli[i];
            let bound = u64::MAX - (u64::MAX % q);
            for j in 0..pout.coeffs[i].len() {
                loop {
                    self.xof.read(&mut buff);
                    let x = u64::from_be_bytes(buff);
                    if x < bound {
                        pout.coeffs[i][j] = x % q;
                        break;
                    }
                }
            }
        }
    }

    /// Reads a monomial.
    /// Output polynomial is in NTT domain.
    pub fn read_monomial(&mut self, ring: &Ring) -> Poly {
        let mut mout = ring.new_poly();
        self.read_monomial_assign(ring, &mut mout);
        return mout;
    }

    /// Reads a monomial.
    /// Output polynomial is in NTT domain.
    pub fn read_monomial_assign(&mut self, ring: &Ring, mout: &mut Poly) {
        mout.clear();

        let mut buff = [0u8; 8];
        let bound = usize::MAX - (usize::MAX % (2 * ring.degree));
        loop {
            self.xof.read(&mut buff);
            let mut x = usize::from_be_bytes(buff);

            if x < bound {
                x %= 2 * ring.degree;

                let s = x & 1;
                if s == 0 {
                    for i in 0..mout.coeffs.len() {
                        mout.coeffs[i][x >> 1] = 1
                    }
                } else {
                    for i in 0..mout.coeffs.len() {
                        mout.coeffs[i][x >> 1] = ring.moduli[i] - 1
                    }
                }

                mout.is_ntt = false;
                ring.ntt(mout);
                return;
            }
        }
    }
}
