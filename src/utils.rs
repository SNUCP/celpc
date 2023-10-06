use primitive_types::U256;

pub fn mod_exp(x: U256, y: usize, m: U256) -> U256 {
    let mut res = U256::from(1);
    let mut x = x;
    let mut y = y;
    while y > 0 {
        if y & 1 == 1 {
            res = (res * x) % m;
        }
        y >>= 1;
        x = (x * x) % m;
    }
    res
}
