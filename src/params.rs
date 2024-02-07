use ethnum::U256;
use primitive_types::U512;
use rug::Integer;

use super::utils::*;

use super::ring::*;

pub const LAMBDA: usize = 128;

#[derive(Clone)]
pub struct Parameters {
    /// mu is the size of commitment key A0, and commitment itself.
    pub mu: usize,
    /// nu is the size of commitment key A1.
    pub nu: usize,
    /// munu equals mu + nu.
    pub munu: usize,

    /// n is the dimension of the message vector.
    /// Equals to l*s.
    pub n: usize,
    /// d is the degree of the ring R_q.
    pub d: usize,
    /// q is the modulus of the ring R_q.
    pub q: u64,

    /// b is the base of the message.
    pub b: u64,
    /// s is the number of "slots".
    /// Essentially, it represents how many Z_p message
    /// is encoded in a single ring element.
    pub s: usize,
    /// l is the length of the packed message.
    /// We often pack n Z_p messages into l ring elements.
    pub l: usize,
    /// kap is number of digits of one message.
    /// Equals to d/s.
    pub kap: usize,
    /// p is the modulus for Z_p.
    /// Equals to b^kap + 1.
    pub p: U256,
    /// p512 is a U512 representation of p.
    pub p512: U512,
    /// pbig is a rug::Integer representation of p.
    pub pbig: Integer,
    /// pr is a barret constant for p.
    pub pr: BarrettConstant,

    /// m is the length of encoded polynomials.
    /// We encode N = m*n degree of polynomial in Z_p,
    /// into m*l number of polynomials in R_q.
    pub m: usize,

    /// ringq is the ring R_q.
    pub ringq: Ring,

    pub s1: f64,
    pub s2: f64,
    pub s3: f64,
    pub sig1: f64,
    pub sig2: f64,
    pub sig3: f64,
}

impl Parameters {
    /// Create a new Parameter.
    pub fn new(
        mu: usize,
        nu: usize,
        d: usize,
        q: u64,
        b: u64,
        s: usize,
        l: usize,
        m: usize,
        s1: f64,
        s2: f64,
        s3: f64,
        sig1: f64,
        sig2: f64,
        sig3: f64,
    ) -> Parameters {
        let kap = d / s;
        let p = U256::from(u128::pow(b as u128, kap as u32) + 1);

        let bbig = Integer::from(b);
        let mut pbig = Integer::from(1);
        for _ in 0..kap {
            pbig *= &bbig;
        }
        pbig += 1;

        Self {
            mu: mu,
            nu: nu,
            munu: mu + nu,
            n: s * l,
            d: d,
            q: q,

            b: b,
            s: s,
            l: l,
            kap: kap,
            p: p,
            p512: U512::from_big_endian(&p.to_be_bytes()),
            pbig: pbig,
            pr: BarrettConstant::new(p),

            m: m,

            ringq: Ring::new(d, &[q]),

            s1: s1,
            s2: s2,
            s3: s3,

            sig1: sig1,
            sig2: sig2,
            sig3: sig3,
        }
    }

    pub fn small() -> Parameters {
        Parameters::new(
            1,
            1,
            1 << 11,
            0xfffffffffffc001,
            248,
            1 << 7,
            1 << 7,
            4,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        )
    }
}
