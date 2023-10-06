use primitive_types::U256;

use super::ring::*;

#[derive(Clone)]
pub struct Parameters {
    /// mu is the size of commitment key A0, and commitment itself.
    pub mu: usize,
    /// nu is the size of commitment key A1.
    pub nu: usize,
    /// munu equals mu + nu.
    pub munu: usize,

    /// n is the size of the ring.
    pub n: usize,
    /// q is the modulus of the ring.
    pub q: u64,
    /// ringq is the ring with modulus q.
    pub ringq: Ring,
    /// ringmul is the ring with modulus qmul.
    pub ringmul: Ring,

    /// b is the base of the message.
    pub b: u64,
    /// m is the number of message
    /// packed into one ring element.
    pub m: usize,
    /// kap is number of digits of one message.
    /// Equals to n/m.
    pub kap: usize,
    /// p is the modulus of the message.
    /// Equals to b^(n/m)+1.
    pub p: U256,

    /// l is the length of packed message.
    pub l: usize,

    /// k is the number of commitments.
    pub k: usize,

    /// r is the number of repeats in Ajtai Protocol.
    pub r: usize,

    /// M is the dimension of R1CS matrix.
    /// Equals to klm.
    pub M: usize,

    pub s1: f64,
    pub s2: f64,
    pub sig1: f64,
    pub sig2: f64,
}

impl Parameters {
    /// Create a new Parameter.
    pub fn new(
        mu: usize,
        nu: usize,
        n: usize,
        q: u64,
        qmul: Vec<u64>,
        b: u64,
        m: usize,
        l: usize,
        k: usize,
        r: usize,
        s1: f64,
        s2: f64,
        sig1: f64,
        sig2: f64,
    ) -> Self {
        Self {
            mu,
            nu,
            munu: mu + nu,
            n,
            q,
            ringq: Ring::new(n, &vec![q]),
            ringmul: Ring::new(4 * k * l * m, &qmul),
            b,
            m,
            kap: n / m,
            p: U256::from(u128::pow(b as u128, (n / m) as u32) + 1),
            l,
            k,
            r,
            M: k * l * m,
            s1,
            s2,
            sig1,
            sig2,
        }
    }

    pub fn small() -> Parameters {
        Parameters::new(
            4,
            4,
            1 << 4,
            0xfffffffffffc001,
            vec![
                0xfffffffa5000001,
                // 0xfffffffaf800001,
                // 0xfffffffba000001,
                // 0xfffffffd9000001,
                // 0xfffffffe0800001,
                // 0xfffffffe8000001,
                // 0xffffffffd000001,
            ],
            14,
            1 << 3,
            1 << 3,
            1 << 3,
            3,
            f64::exp(2.98),
            f64::exp(12.73),
            f64::exp(3.98),
            f64::exp(5.77),
        )
    }
}

impl Default for Parameters {
    fn default() -> Self {
        Self::new(
            1,
            1,
            1 << 11,
            0xfffffffffffc001,
            vec![
                0xfffffffa5000001,
                0xfffffffaf800001,
                0xfffffffba000001,
                0xfffffffd9000001,
                0xfffffffe0800001,
                0xfffffffe8000001,
                0xffffffffd000001,
            ],
            248,
            1 << 7,
            1 << 6,
            1 << 7,
            11,
            f64::exp(2.98),
            f64::exp(12.73),
            f64::exp(3.98),
            f64::exp(5.77),
        )
    }
}
