use ethnum::U256;

use crate::utils::*;

pub struct SparseMatrix {
    pub n: usize,
    pub q: U256,
    pub qbmod: BarrettConstant,

    // CSR Format
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<usize>,
    pub values: Vec<U256>,
}

impl SparseMatrix {
    /// Create a new SparseMatrix with n rows and q modulus.
    pub fn new(n: usize, q: U256) -> SparseMatrix {
        SparseMatrix {
            n,
            q: q,
            qbmod: BarrettConstant::new(q),

            row_ptr: vec![0; n + 1],
            col_idx: vec![],
            values: vec![],
        }
    }

    /// Create a identity Matrix with n rows and q modulus.
    pub fn new_identity(n: usize, q: U256) -> SparseMatrix {
        let mut row_ptr = vec![0; n + 1];
        let mut col_idx = vec![0; n];
        let mut values = vec![U256::ZERO; n];

        for i in 0..n {
            row_ptr[i] = i;
            col_idx[i] = i;
            values[i] = U256::ONE;
        }
        row_ptr[n] = n;

        SparseMatrix {
            n: n,
            q: q,
            qbmod: BarrettConstant::new(q),

            row_ptr: row_ptr,
            col_idx: col_idx,
            values: values,
        }
    }

    /// Transforms a sparse matrix to dense matrix.
    pub fn to_dense(&self) -> Vec<Vec<U256>> {
        let mut m = vec![vec![U256::ZERO; self.n]; self.n];

        for i in 0..self.n {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                m[i][self.col_idx[j]] = self.values[j];
            }
        }

        return m;
    }

    /// Multiplies v and returns the result.
    pub fn mul_vec(&self, v: &[U256]) -> Vec<U256> {
        let mut vout = vec![U256::ZERO; self.n];
        self.mul_vec_assign(v, &mut vout);
        return vout;
    }

    /// Multiplies v and writes it to vout.
    pub fn mul_vec_assign(&self, v: &[U256], vout: &mut [U256]) {
        for i in 0..self.n {
            vout[i] = U256::ZERO;
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                vout[i] += bmod(self.values[j] * v[self.col_idx[j]], self.q, self.qbmod);
                if vout[i] >= self.q {
                    vout[i] -= self.q;
                }
            }
        }
    }

    /// Multiplies v and adds it to vout.
    pub fn mul_vec_add_assign(&self, v: &[U256], vout: &mut [U256]) {
        for i in 0..self.n {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                vout[i] += bmod(self.values[j] * v[self.col_idx[j]], self.q, self.qbmod);
                if vout[i] >= self.q {
                    vout[i] -= self.q;
                }
            }
        }
    }

    /// Multiplies v and subtracts it from vout.
    pub fn mul_vec_sub_assign(&self, v: &[U256], vout: &mut [U256]) {
        for i in 0..self.n {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                let tmp = bmod(self.values[j] * v[self.col_idx[j]], self.q, self.qbmod);
                if vout[i] >= tmp {
                    vout[i] = vout[i] - tmp
                } else {
                    vout[i] = vout[i] + self.q - tmp
                }
            }
        }
    }

    /// Transposes and multiplies v and returns the result.
    pub fn transpose_mul_vec(&self, v: &[U256]) -> Vec<U256> {
        let mut vout = vec![U256::ZERO; self.n];
        self.transpose_mul_vec_assign(v, &mut vout);
        return vout;
    }

    /// Transposes and multiplies v and writes it to vout.
    pub fn transpose_mul_vec_assign(&self, v: &[U256], vout: &mut [U256]) {
        for i in 0..self.n {
            vout[i] = U256::ZERO;
        }

        for i in 0..self.n {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                vout[self.col_idx[j]] += bmod(self.values[j] * v[i], self.q, self.qbmod);
                if vout[self.col_idx[j]] >= self.q {
                    vout[self.col_idx[j]] -= self.q;
                }
            }
        }
    }

    /// Transposes and multiplies v and adds it to vout.
    pub fn transpose_mul_vec_add_assign(&self, v: &[U256], vout: &mut [U256]) {
        for i in 0..self.n {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                vout[self.col_idx[j]] += bmod(self.values[j] * v[i], self.q, self.qbmod);
                if vout[self.col_idx[j]] >= self.q {
                    vout[self.col_idx[j]] -= self.q;
                }
            }
        }
    }

    /// Transposes and multiplies v and subtracts it from vout.
    pub fn transpose_mul_vec_sub_assign(&self, v: &[U256], vout: &mut [U256]) {
        for i in 0..self.n {
            for j in self.row_ptr[i]..self.row_ptr[i + 1] {
                let tmp = bmod(self.values[j] * v[i], self.q, self.qbmod);
                if vout[self.col_idx[j]] >= tmp {
                    vout[self.col_idx[j]] = vout[self.col_idx[j]] - tmp
                } else {
                    vout[self.col_idx[j]] = vout[self.col_idx[j]] + self.q - tmp;
                }
            }
        }
    }
}
