use ethnum::U256;
use primitive_types::U512;

#[inline]
pub fn mod_exp(x: U256, y: usize, m: U512) -> U256 {
    let mut res = U512::from(1);
    let mut x = mod_up(x);
    let mut y = y;
    while y > 0 {
        if y & 1 == 1 {
            res = (res * x) % m;
        }
        y >>= 1;
        x = (x * x) % m;
    }
    return mod_down(res);
}

#[inline]
pub fn mod_up(x: U256) -> U512 {
    return U512::from_big_endian(&x.to_be_bytes());
}

#[inline]
pub fn mod_down(x: U512) -> U256 {
    let mut res = U256::ZERO;
    res.0[0] = ((x.0[1] as u128) << 64) | (x.0[0] as u128);
    res.0[1] = ((x.0[3] as u128) << 64) | (x.0[2] as u128);
    return res;
}

#[inline]
pub fn inner_product(x: &[U256], y: &[U256], q: U512) -> U256 {
    let mut res = U512::zero();
    for i in 0..x.len() {
        res += mod_up(x[i]) * mod_up(y[i]);
        res %= q;
    }
    return mod_down(res);
}
