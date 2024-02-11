#[derive(Debug, Clone)]
pub struct Poly {
    pub coeffs: Vec<Vec<u64>>,
    pub is_ntt: bool,
}

impl Poly {
    /// Clears the polynomial.
    /// Does not change the NTT flag.
    pub fn clear(&mut self) {
        for i in 0..self.coeffs.len() {
            for j in 0..self.coeffs[i].len() {
                self.coeffs[i][j] = 0;
            }
        }
    }

    /// Checks if two polynomials are equal.
    pub fn equal(&self, other: &Poly) -> bool {
        if self.coeffs.len() != other.coeffs.len() {
            return false;
        }

        for i in 0..self.coeffs.len() {
            if self.coeffs[i].len() != other.coeffs[i].len() {
                return false;
            }
            for j in 0..self.coeffs[i].len() {
                if self.coeffs[i][j] != other.coeffs[i][j] {
                    return false;
                }
            }
        }

        return true;
    }

    /// is_zero checks if the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        for i in 0..self.coeffs.len() {
            for j in 0..self.coeffs[i].len() {
                if self.coeffs[i][j] != 0 {
                    return false;
                }
            }
        }
        return true;
    }
}
