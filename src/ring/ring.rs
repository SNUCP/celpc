use super::Poly;
use rug::{Complete, Integer};
use tfhe_ntt::prime64::Plan;

/// Ring is a polynomial ring Z_q[X]/(X^N + 1),
/// where N is a power of two.
#[derive(Clone)]
pub struct Ring {
    degree: usize,
    modulus: Vec<u64>,
    plan: Vec<Plan>,

    modulus_big: Integer,
    rns_gadget: Vec<Integer>,
}

impl Ring {
    /// Creates a ring.
    pub fn new(n: usize, q: &[u64]) -> Self {
        let plan = q.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();

        let q_full = q.iter().fold(Integer::from(1), |acc, &q| acc * q);
        let mut rns_gadget = vec![Integer::ZERO; q.len()];
        for i in 0..q.len() {
            let qi = Integer::from(q[i]);
            let q_div = (&q_full / &qi).complete();
            let q_div_inv = q_div.clone().invert(&qi).unwrap();
            rns_gadget[i] = q_div * q_div_inv;
        }

        Self {
            degree: n,
            modulus: q.to_vec(),
            plan,

            modulus_big: q_full,
            rns_gadget: rns_gadget,
        }
    }

    /// degree returns the degree of the ring.
    pub fn degree(&self) -> usize {
        return self.degree;
    }

    /// level returns the number of moduli in the ring.
    pub fn level(&self) -> usize {
        return self.modulus.len();
    }

    /// modulus returns the modulus of the ring.
    pub fn modulus(&self) -> Vec<u64> {
        return self.modulus.clone();
    }

    /// modulus_idx returns the i-th modulus of the ring.
    pub fn modulus_idx(&self, i: usize) -> u64 {
        return self.modulus[i];
    }

    /// modulus_big returns the modulus of the ring as a big integer.
    pub fn modulus_big(&self) -> Integer {
        return self.modulus_big.clone();
    }

    /// new_poly creates a new polynomial.
    pub fn new_poly(&self) -> Poly {
        return Poly::new(self.degree, self.modulus.len());
    }

    /// ntt_inplace computes the ntt of p in-place.
    pub fn ntt_inplace(&self, p: &mut Poly) {
        for (i, plan) in self.plan.iter().enumerate() {
            plan.fwd(&mut p.coeffs[i]);
        }
    }

    /// intt_inplace computes the intt of p in-place.
    pub fn intt_inplace(&self, p: &mut Poly) {
        for (i, plan) in self.plan.iter().enumerate() {
            plan.inv(&mut p.coeffs[i]);
            plan.normalize(&mut p.coeffs[i]);
        }
    }

    /// add returns p_out = p0 + p1.
    pub fn add(&self, p0: &Poly, p1: &Poly) -> Poly {
        let mut p_out = Poly::new(self.degree, self.modulus.len());
        self.add_assign(p0, p1, &mut p_out);
        return p_out;
    }

    /// add_assign computes p_out = p0 + p1.
    pub fn add_assign(&self, p0: &Poly, p1: &Poly, p_out: &mut Poly) {
        for i in 0..self.modulus.len() {
            let qi = self.modulus[i];
            for j in 0..self.degree {
                p_out.coeffs[i][j] = p0.coeffs[i][j] + p1.coeffs[i][j];
                if p_out.coeffs[i][j] >= qi {
                    p_out.coeffs[i][j] -= qi;
                }
            }
        }
    }

    /// add_inplace computes p_out += p0.
    pub fn add_inplace(&self, p0: &Poly, p_out: &mut Poly) {
        for i in 0..self.modulus.len() {
            let qi = self.modulus[i];
            for j in 0..self.degree {
                p_out.coeffs[i][j] += p0.coeffs[i][j];
                if p_out.coeffs[i][j] >= qi {
                    p_out.coeffs[i][j] -= qi;
                }
            }
        }
    }

    /// sub returns p_out = p0 - p1.
    pub fn sub(&self, p0: &Poly, p1: &Poly) -> Poly {
        let mut p_out = Poly::new(self.degree, self.modulus.len());
        self.sub_assign(p0, p1, &mut p_out);
        return p_out;
    }

    /// sub_assign computes p_out = p0 - p1.
    pub fn sub_assign(&self, p0: &Poly, p1: &Poly, p_out: &mut Poly) {
        for i in 0..self.modulus.len() {
            let qi = self.modulus[i];
            for j in 0..self.degree {
                p_out.coeffs[i][j] = p0.coeffs[i][j].wrapping_sub(p1.coeffs[i][j]);
                if p_out.coeffs[i][j] >= qi {
                    p_out.coeffs[i][j] += qi;
                }
            }
        }
    }

    /// sub_inplace computes p_out -= p0.
    pub fn sub_inplace(&self, p0: &Poly, p_out: &mut Poly) {
        for i in 0..self.modulus.len() {
            let qi = self.modulus[i];
            for j in 0..self.degree {
                p_out.coeffs[i][j] = p_out.coeffs[i][j].wrapping_sub(p0.coeffs[i][j]);
                if p_out.coeffs[i][j] >= qi {
                    p_out.coeffs[i][j] -= qi;
                }
            }
        }
    }

    /// mul returns p_out = p0 * p1.
    pub fn mul(&self, p0: &Poly, p1: &Poly) -> Poly {
        let mut p_out = Poly::new(self.degree, self.modulus.len());
        self.mul_assign(p0, p1, &mut p_out);
        return p_out;
    }

    /// mul_assign computes p_out = p0 * p1.
    pub fn mul_assign(&self, p0: &Poly, p1: &Poly, p_out: &mut Poly) {
        for (i, plan) in self.plan.iter().enumerate() {
            p_out.coeffs[i].fill(0);
            plan.mul_accumulate(&mut p_out.coeffs[i], &p0.coeffs[i], &p1.coeffs[i]);
        }
    }

    /// mul_add_assign computes p_out += p0 * p1.
    pub fn mul_add_assign(&self, p0: &Poly, p1: &Poly, p_out: &mut Poly) {
        for (i, plan) in self.plan.iter().enumerate() {
            plan.mul_accumulate(&mut p_out.coeffs[i], &p0.coeffs[i], &p1.coeffs[i]);
        }
    }

    /// mul_sub_assign computes p_out -= p0 * p1.
    pub fn mul_sub_assign(&self, p0: &Poly, p1: &Poly, p_out: &mut Poly) {
        for (i, plan) in self.plan.iter().enumerate() {
            let qi = self.modulus[i];
            for j in 0..self.degree {
                p_out.coeffs[i][j] = qi - p_out.coeffs[i][j];
            }
            plan.mul_accumulate(&mut p_out.coeffs[i], &p0.coeffs[i], &p1.coeffs[i]);
            for j in 0..self.degree {
                p_out.coeffs[i][j] = qi - p_out.coeffs[i][j];
            }
        }
    }

    fn to_balanced(&self, x: u64, q: u64) -> i64 {
        if x >= q / 2 {
            return x as i64 - q as i64;
        }
        return x as i64;
    }

    /// to_bigint_balanced converts p to a big integer in the range [-q/2, q/2).
    pub fn to_bigint_balanced(&self, p: &Poly) -> Vec<Integer> {
        let mut p_intt = p.clone();
        self.intt_inplace(&mut p_intt);

        let mut big_out = Vec::with_capacity(self.degree);

        for i in 0..self.degree {
            let mut is_small = true;
            let c_i64 = self.to_balanced(p_intt.coeffs[0][i], self.modulus[0]);
            for j in 1..self.modulus.len() {
                let c_j64 = self.to_balanced(p_intt.coeffs[j][i], self.modulus[j]);
                if c_i64 != c_j64 {
                    is_small = false;
                    break;
                }
            }

            if is_small {
                big_out.push(Integer::from(c_i64));
                continue;
            }

            big_out.push(Integer::ZERO);
            for j in 0..self.modulus.len() {
                let mut c = Integer::from(p_intt.coeffs[j][i]);
                c *= &self.rns_gadget[j];
                big_out[i] += &c;
            }
            big_out[i].modulo_mut(&self.modulus_big);
        }

        return big_out;
    }

    /// norm returns the square of 2-norm of p.
    pub fn norm_sq(&self, p: &Poly) -> f64 {
        let big = self.to_bigint_balanced(p);

        let mut norm = Integer::ZERO;
        for i in 0..self.degree {
            norm += &big[i] * &big[i];
        }

        return norm.to_f64();
    }
}
