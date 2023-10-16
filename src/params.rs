use ethnum::U256;

use crate::utils::*;

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
    /// pbmod is barret constant for
    /// modulo p.
    pub pbmod: BarrettConstant,

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
    ) -> Parameters {
        let kap = n / m;
        let p = U256::from(u128::pow(b as u128, kap as u32) + 1);
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
            kap: kap,
            p: p,
            pbmod: BarrettConstant::new(p),
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

    pub fn bit_18() -> Parameters {
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
            ],
            248,
            1 << 7, // m
            1 << 5, // l
            1 << 6, // k
            9,      // r
            f64::exp(2.98),
            f64::exp(12.73),
            f64::exp(3.98),
            f64::exp(5.77),
        )
    }

    pub fn bit_20() -> Parameters {
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
            ],
            248,
            1 << 7, // m
            1 << 6, // l
            1 << 7, // k
            9,      // r
            f64::exp(2.98),
            f64::exp(12.73),
            f64::exp(3.98),
            f64::exp(5.77),
        )
    }

    pub fn bit_22() -> Parameters {
        Self::new(
            1,
            1,
            1 << 11,
            0xfffffffffffc001,
            vec![
                0xffffffee8000001,
                0xfffffff16000001,
                0xfffffff4e000001,
                0xfffffffba000001,
                0xfffffffe8000001,
            ],
            248,
            1 << 7, // m
            1 << 7, // l
            1 << 8, // k
            9,      // r
            f64::exp(2.98),
            f64::exp(12.73),
            f64::exp(3.98),
            f64::exp(5.77),
        )
    }
}
