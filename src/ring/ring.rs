use super::*;
use concrete_ntt::prime64::Plan;
use ethnum::U256;
use rug::ops::*;
use rug::Integer;

#[derive(Clone)]
pub struct Ring {
    pub degree: usize,

    pub moduli: Vec<u64>,
    pub moduli_big: Integer,
    pub gadget: Vec<Integer>,

    pub plans: Vec<Plan>,
}

impl Ring {
    /// Create a new ring with the given degree and moduli.
    pub fn new(degree: usize, moduli: &[u64]) -> Ring {
        let mut plans = Vec::with_capacity(moduli.len());
        for &q in moduli {
            let plan = Plan::try_new(degree, q).unwrap();
            plans.push(plan);
        }

        let mut gadget = vec![Integer::from(0); degree];
        let qfull = moduli.iter().fold(Integer::from(1), |acc, x| acc * x);
        for (i, &q) in moduli.iter().enumerate() {
            let qbig = Integer::from(q);

            let qstar = qfull.clone().div_exact(&qbig);
            let qinv = qstar.clone().invert(&qbig).unwrap();
            gadget[i] = qstar * qinv;
        }

        return Ring {
            degree: degree,

            moduli: moduli.to_vec(),
            moduli_big: qfull,
            gadget: gadget,

            plans: plans,
        };
    }

    /// Create a new polynomial with the given degree and modulus.
    /// The polynomial is not in NTT form.
    pub fn new_poly(&self) -> Poly {
        let coeffs = vec![vec![0; self.degree]; self.moduli.len()];
        Poly {
            coeffs: coeffs,
            is_ntt: false,
        }
    }

    /// Create a new polynomial with the given degree and modulus.
    /// The polynomial is in NTT form.
    pub fn new_ntt_poly(&self) -> Poly {
        let coeffs = vec![vec![0; self.degree]; self.moduli.len()];
        Poly {
            coeffs: coeffs,
            is_ntt: true,
        }
    }

    /// Applies NTT to the given polynomial.
    /// Does nothing if the polynomial is already in NTT form.
    pub fn ntt(&self, poly: &mut Poly) {
        if poly.is_ntt {
            panic!("Polynomial is already in NTT form.");
        }

        for i in 0..self.moduli.len() {
            self.plans[i].fwd(&mut poly.coeffs[i]);
        }
        poly.is_ntt = true;
    }

    /// Applies Inverse NTT to the given polynomial.
    /// Does nothing if the polynomial is not in NTT form.
    pub fn intt(&self, poly: &mut Poly) {
        if !poly.is_ntt {
            panic!("Polynomial is not in NTT form.");
        }

        for i in 0..self.moduli.len() {
            self.plans[i].inv(&mut poly.coeffs[i]);
            self.plans[i].normalize(&mut poly.coeffs[i]);
        }
        poly.is_ntt = false;
    }

    /// Sets the coefficient of the polynomial at the given index to the given value.
    pub fn set_coeff(&self, p: &mut Poly, idx: usize, coeff: U256) {
        for (i, &q) in self.moduli.iter().enumerate() {
            p.coeffs[i][idx] = (coeff % U256::from(q)).as_u64();
        }
    }

    /// Gets the coefficient of the polynomial at the given index as Integer.
    pub fn get_coeff(&self, p: &Poly, idx: usize) -> Integer {
        let mut coeff = Integer::from(0);
        for i in 0..self.moduli.len() {
            coeff += Integer::from(p.coeffs[i][idx]) * &self.gadget[i];
        }
        coeff.rem_euc_assign(&self.moduli_big);
        return coeff;
    }

    /// Gets the balanced coefficient of the polynomial at the given index as Integer.
    pub fn get_coeff_balanced(&self, p: &Poly, idx: usize) -> Integer {
        let coeff = self.get_coeff(p, idx);
        let qhalf = self.moduli_big.clone() >> 1;
        if coeff >= qhalf {
            coeff - &self.moduli_big
        } else {
            coeff
        }
    }

    /// Maps a polynomial in R_Q to balanced form.
    pub fn to_balanced(&self, p: &Poly) -> Vec<Vec<i64>> {
        let mut vout = vec![vec![0; self.degree]; self.moduli.len()];
        self.to_balanced_assign(p, &mut vout);
        return vout;
    }

    /// Maps a polynomial in R_Q to balanced form.
    pub fn to_balanced_assign(&self, p: &Poly, vout: &mut Vec<Vec<i64>>) {
        let mut buff = p.clone();
        if p.is_ntt {
            self.intt(&mut buff);
        }

        for i in 0..self.moduli.len() {
            let qs = self.moduli[i] as i64;
            for (j, &c) in buff.coeffs[i].iter().enumerate() {
                vout[i][j] = c as i64;
                if vout[i][j] > qs / 2 {
                    vout[i][j] -= qs;
                }
            }
        }
    }

    /// Adds p1, p2 and returns it.
    pub fn add(&self, p1: &Poly, p2: &Poly) -> Poly {
        if p1.is_ntt != p2.is_ntt {
            panic!("Polynomials must be in the same form.");
        }

        let mut pout = self.new_poly();
        self.add_assign(p1, p2, &mut pout);
        return pout;
    }

    /// Adds p1, p2 and writes it to pout.
    pub fn add_assign(&self, p1: &Poly, p2: &Poly, pout: &mut Poly) {
        if p1.is_ntt != p2.is_ntt {
            panic!("Polynomials must be in the same form.");
        }

        pout.is_ntt = p1.is_ntt;
        for (i, &q) in self.moduli.iter().enumerate() {
            for j in 0..self.degree {
                pout.coeffs[i][j] = p1.coeffs[i][j] + p2.coeffs[i][j];
                if pout.coeffs[i][j] >= q {
                    pout.coeffs[i][j] -= q;
                }
            }
        }
    }

    /// Adds p to pout.
    pub fn add_inplace(&self, p: &Poly, pout: &mut Poly) {
        if p.is_ntt != pout.is_ntt {
            panic!("Polynomials must be in the same form.");
        }

        for (i, &q) in self.moduli.iter().enumerate() {
            for j in 0..self.degree {
                pout.coeffs[i][j] += p.coeffs[i][j];
                if pout.coeffs[i][j] >= q {
                    pout.coeffs[i][j] -= q;
                }
            }
        }
    }

    // Subtracts p2 from p1 and returns it.
    pub fn sub(&self, p1: &Poly, p2: &Poly) -> Poly {
        if p1.is_ntt != p2.is_ntt {
            panic!("Polynomials must be in the same form.");
        }

        let mut pout = self.new_poly();
        self.sub_assign(p1, p2, &mut pout);
        return pout;
    }

    // Subtracts p2 from p1 and writes it to pout.
    pub fn sub_assign(&self, p1: &Poly, p2: &Poly, pout: &mut Poly) {
        if p1.is_ntt != p2.is_ntt {
            panic!("Polynomials must be in the same form.");
        }

        pout.is_ntt = p1.is_ntt;
        for (i, &q) in self.moduli.iter().enumerate() {
            for j in 0..self.degree {
                if p1.coeffs[i][j] >= p2.coeffs[i][j] {
                    pout.coeffs[i][j] = p1.coeffs[i][j] - p2.coeffs[i][j];
                } else {
                    pout.coeffs[i][j] = p1.coeffs[i][j] + (q - p2.coeffs[i][j]);
                }
            }
        }
    }

    /// Subtracts p from pout.
    pub fn sub_inplace(&self, p: &Poly, pout: &mut Poly) {
        if p.is_ntt != pout.is_ntt {
            panic!("Polynomials must be in the same form.");
        }

        for (i, &q) in self.moduli.iter().enumerate() {
            for j in 0..self.degree {
                if pout.coeffs[i][j] > p.coeffs[i][j] {
                    pout.coeffs[i][j] -= p.coeffs[i][j];
                } else {
                    pout.coeffs[i][j] += q - p.coeffs[i][j];
                }
            }
        }
    }

    /// Multiplies p1, p2 and returns it.
    /// The input polynomials are assumed to be in NTT form.
    /// The output polynomial is in NTT form.
    pub fn mul(&self, p1: &Poly, p2: &Poly) -> Poly {
        let mut pout = self.new_ntt_poly();
        self.mul_assign(p1, p2, &mut pout);
        return pout;
    }

    /// Multiplies p1, p2 and writes it to pout.
    /// The input polynomials are assumed to be in NTT form.
    /// The output polynomial is in NTT form.
    pub fn mul_assign(&self, p1: &Poly, p2: &Poly, pout: &mut Poly) {
        if !(p1.is_ntt && p2.is_ntt && pout.is_ntt) {
            panic!("Polynomials must be in NTT form.");
        }

        for i in 0..self.moduli.len() {
            pout.coeffs[i].fill(0);
            self.plans[i].mul_accumulate(&mut pout.coeffs[i], &p1.coeffs[i], &p2.coeffs[i]);
        }
    }

    /// Multiplies p to pout.
    /// The input polynomials are assumed to be in NTT form.
    /// The output polynomial is in NTT form.
    pub fn mul_inplace(&self, p: &Poly, pout: &mut Poly) {
        if !(p.is_ntt && pout.is_ntt) {
            panic!("Polynomials must be in NTT form.");
        }

        let mut buff = vec![0u64; self.degree];
        for i in 0..self.moduli.len() {
            buff.fill(0);
            self.plans[i].mul_accumulate(&mut buff, &p.coeffs[i], &pout.coeffs[i]);
            buff.clone_into(&mut pout.coeffs[i]);
        }
    }

    /// Multiplies p1, p2 and adds it to pout.
    /// The input polynomials are assumed to be in NTT form.
    /// The output polynomial is in NTT form.
    pub fn mul_add_inplace(&self, p1: &Poly, p2: &Poly, pout: &mut Poly) {
        if !(p1.is_ntt && p2.is_ntt && pout.is_ntt) {
            panic!("Polynomials must be in NTT form.");
        }

        for i in 0..self.moduli.len() {
            self.plans[i].mul_accumulate(&mut pout.coeffs[i], &p1.coeffs[i], &p2.coeffs[i]);
        }
    }

    /// Multiplies p1, p2 and subracts it from pout.
    /// The input polynomials are assumed to be in NTT form.
    /// The output polynomial is in NTT form.
    pub fn mul_sub_inplace(&self, p1: &Poly, p2: &Poly, pout: &mut Poly) {
        if !(p1.is_ntt && p2.is_ntt && pout.is_ntt) {
            panic!("Polynomials must be in NTT form.");
        }

        let mut buff = vec![0u64; self.degree];
        for i in 0..self.moduli.len() {
            for j in 0..self.degree {
                buff[j] = self.moduli[i] - p1.coeffs[i][j];
            }
            self.plans[i].mul_accumulate(&mut pout.coeffs[i], &buff, &p2.coeffs[i]);
        }
    }

    /// Returns the pair of (L2 norm^2, L-infinity norm) of the given polynomial.
    pub fn norm(&self, p: &Poly) -> (f64, f64) {
        let buff = self.to_balanced(p);
        let mut l2 = 0.0;
        let mut linf = 0.0;

        for i in 0..buff.len() {
            for j in 0..buff[i].len() {
                let c = (buff[i][j] as f64).abs();
                l2 += c * c;
                if c > linf {
                    linf = c;
                }
            }
        }

        return (l2, linf);
    }
}
