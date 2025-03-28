/// Poly is a RNS polynomial.
#[derive(Clone)]
pub struct Poly {
    pub coeffs: Vec<Vec<u64>>,
}

impl Poly {
    /// Creates a polynomial.
    pub fn new(n: usize, level: usize) -> Poly {
        Poly {
            coeffs: vec![vec![0; n]; level],
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

    /// Clears the polynomial.
    pub fn clear(&mut self) {
        for i in 0..self.coeffs.len() {
            self.coeffs[i].fill(0);
        }
    }
}
