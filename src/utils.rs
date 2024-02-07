use ethnum::U256;
use primitive_types::U512;

#[inline]
pub fn mod_exp(x: U256, y: usize, m: U256, r: BarrettConstant) -> U256 {
    let mut res = U256::ONE;
    let mut x = x;
    let mut y = y;
    while y > 0 {
        if y & 1 == 1 {
            res = bmod(res * x, m, r);
        }
        y >>= 1;
        x = bmod(x * x, m, r);
    }
    res
}

#[derive(Clone, Copy, Debug)]
pub struct BarrettConstant(U256, U256);

impl BarrettConstant {
    pub fn new(q: U256) -> BarrettConstant {
        let r = U256::MAX / q;
        let r0 = U256([r.0[0], 0]);
        let r1 = U256([r.0[1], 0]);

        BarrettConstant(r0, r1)
    }
}

#[inline]
pub fn mulmod(x: U256, y: U256, q: U512) -> U256 {
    let x512 = U512::from_big_endian(&x.to_be_bytes());
    let y512 = U512::from_big_endian(&y.to_be_bytes());
    let xy = x512 * y512 % q;

    let mut res = U256::ZERO;
    res.0[0] = ((xy.0[1] as u128) << 64) | (xy.0[0] as u128);
    res.0[1] = ((xy.0[3] as u128) << 64) | (xy.0[2] as u128);
    return res;
}

#[inline]
pub fn bmod(x: U256, q: U256, r: BarrettConstant) -> U256 {
    let x0 = U256([x.0[0], 0]);
    let x1 = U256([x.0[1], 0]);

    let t = r.1 * x1 + ((r.1 * x0 + r.0 * x1 + ((r.0 * x0) >> 128)) >> 128);
    let xr = x - t * q;
    if xr >= q {
        xr - q
    } else {
        xr
    }
}

#[inline]
pub fn inner_product(x: &[U256], y: &[U256], q: U256, r: BarrettConstant) -> U256 {
    let mut res = U256::ZERO;
    for i in 0..x.len() {
        res += bmod(x[i] * y[i], q, r);
    }
    bmod(res, q, r)
}

/// This reverses x!
#[inline]
pub fn inner_product_rev(x: &[U256], y: &[U256], q: U256, r: BarrettConstant) -> U256 {
    let mut res = U256::ZERO;
    for i in 0..x.len() {
        res += bmod(x[x.len() - i - 1] * y[i], q, r);
    }
    bmod(res, q, r)
}

#[inline]
pub fn hadamard_product(x: &[U256], y: &[U256], q: U256, r: BarrettConstant, vout: &mut [U256]) {
    for i in 0..x.len() {
        vout[i] = bmod(x[i] * y[i], q, r);
    }
}

#[inline]
pub fn hadamard_product_inplace(x: &[U256], q: U256, r: BarrettConstant, vout: &mut [U256]) {
    for i in 0..x.len() {
        vout[i] = bmod(x[i] * vout[i], q, r);
    }
}
