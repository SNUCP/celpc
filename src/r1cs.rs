use primitive_types::U256;
use rug::{ops::RemRounding, Integer};

use crate::{csprng::*, ring::*, *};

pub struct R1CSProof {
    pub dcommit: Vec<Vec<Poly>>,
    pub ecommit0: Vec<Poly>,
    pub ecommit1: Vec<Poly>,
    pub ecommit2: Vec<Poly>,
    pub hcommit: Vec<Vec<Poly>>,

    pub t: Vec<Poly>,
    pub tau: Vec<Poly>,
    pub d: Vec<Poly>,
    pub delta: Vec<Poly>,
    pub h: Vec<Poly>,
    pub eta: Vec<Poly>,
}

pub struct R1CSProver<'a> {
    pub params: &'a Parameters,
    pub encoder: Encoder<'a>,
    pub sampler: KarneySampler,
    pub oracle: Oracle,

    pub key: &'a CommitKey,

    pub comitter: Comitter<'a>,
    pub ajprover: AjtaiProver<'a>,
    pub ajverifier: AjtaiVerifier<'a>,

    pub e_msg_stddev1: f64,
    pub e_rand_stddev1: f64,
    pub e_msg_stddev2: f64,
    pub e_rand_stddev2: f64,
}

impl<'a> R1CSProver<'a> {
    /// Create a new R1CSProver.
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> R1CSProver<'a> {
        let stddev1 =
            (params.k as f64).sqrt() * (params.kap as f64) * ((params.b + 1) as f64) * 0.5;
        let stddev2 = ((2 * params.k + 2) as f64).sqrt()
            * (params.kap as f64)
            * ((params.b + 1) as f64)
            * 0.5;

        R1CSProver {
            params: params,
            encoder: Encoder::new(params),
            sampler: KarneySampler::new(),
            oracle: Oracle::new(),

            key: key,

            comitter: Comitter::new(params, key),
            ajprover: AjtaiProver::new(params, key),
            ajverifier: AjtaiVerifier::new(params, key),

            e_msg_stddev1: stddev1 * params.s1,
            e_rand_stddev1: stddev1 * params.sig1,
            e_msg_stddev2: stddev2 * params.s2,
            e_rand_stddev2: stddev2 * params.sig2,
        }
    }

    /// Proves a R1CS instance.
    pub fn prove(
        &mut self,
        A: &SparseMatrix,
        B: &SparseMatrix,
        C: &SparseMatrix,
        a: &[U256],
        b: &[U256],
        c: &[U256],
        t: &[Vec<Poly>],
        tau: &[Vec<Poly>],
        tcommit: &[Vec<Poly>],
    ) -> R1CSProof {
        let params = self.params;

        // Line 2
        // TODO: This can be removed by just accepting decoded t.
        let mut tdec = vec![U256::from(0); params.M];
        for (i, chunk) in tdec.chunks_exact_mut(params.l * params.m).enumerate() {
            self.encoder.decode_chunk_assign(&t[i], chunk)
        }

        // Line 3
        // TODO: Incorrect for now, must feed more information about A, B, C...
        A.values
            .iter()
            .chain(&B.values)
            .chain(&C.values)
            .for_each(|&x| self.oracle.write_u256(x));
        a.iter()
            .chain(b)
            .chain(c)
            .for_each(|&x| self.oracle.write_u256(x));
        tcommit
            .iter()
            .flatten()
            .for_each(|p| self.oracle.write_poly(p));

        // Line 4
        let mut v = vec![U256::from(1); params.M];
        self.oracle.finalize();
        v[1] = self.oracle.read_range(params.p);
        self.oracle.write_u256(v[1]);
        for i in 2..params.M {
            v[i] = (v[i - 1] * v[1]) % params.p;
        }

        // Line 5
        let mut d = vec![U256::from(0); params.M];
        let mut buff_vec = vec![U256::from(0); params.M];

        // B.T*V*a
        for i in 0..params.M {
            buff_vec[i] = (a[i] * v[i]) % params.p; // V*a
        }
        B.transpose_mul_vec_add_assign(&buff_vec, &mut d); // B.T*V*a

        // A.T*V*b
        for i in 0..params.M {
            buff_vec[i] = (b[i] * v[i]) % params.p; // V*b
        }
        A.transpose_mul_vec_add_assign(&buff_vec, &mut d); // B.T*V*a + A.T*V*b
        C.transpose_mul_vec_sub_assign(&v, &mut d); // B.T*V*a + A.T*V*b - C.T*v
        let d_reuse = d.clone(); // Reuse this in g1

        // A.T*V*B*t
        B.mul_vec_assign(&tdec, &mut buff_vec);
        for i in 0..params.M {
            buff_vec[i] = (buff_vec[i] * v[i]) % params.p; // V*B*t
        }
        A.transpose_mul_vec_add_assign(&buff_vec, &mut d); // B.T*V*a + A.T*V*b - C.T*v + A.T*V*B*t

        let mut denc = vec![vec![params.ringq.new_ntt_poly(); params.l]; params.k];
        let mut delta = vec![vec![params.ringq.new_ntt_poly(); params.munu]; params.k];
        let mut dcommit = vec![vec![params.ringq.new_ntt_poly(); params.mu]; params.k];
        for (i, di) in d.chunks_exact(params.l * params.m).enumerate() {
            // Line 7
            self.encoder
                .encode_randomized_chunk_assign(&di, params.s1, &mut denc[i]);

            delta[i].iter_mut().for_each(|p| {
                self.sampler
                    .sample_poly_assign(&params.ringq, 0.0, params.sig1, p)
            });

            // Line 8
            self.comitter
                .commit_assign(&denc[i], &delta[i], &mut dcommit[i]);

            dcommit[i].iter().for_each(|p| self.oracle.write_poly(p));
        }

        // Line 9
        let mut w = vec![U256::from(1); params.M];
        self.oracle.finalize();
        w[1] = self.oracle.read_range(params.p);
        self.oracle.write_u256(w[1]);

        // Line 10
        for i in 2..params.M {
            w[i] = (w[i - 1] * w[1]) % params.p;
        }

        // Line 12
        let mut e0 = vec![params.ringq.new_ntt_poly(); params.l];
        e0.iter_mut().for_each(|p| {
            self.sampler
                .sample_poly_assign(&params.ringq, 0.0, self.e_msg_stddev1, p)
        });
        let mut e1 = vec![params.ringq.new_ntt_poly(); params.l];
        e1.iter_mut().for_each(|p| {
            self.sampler
                .sample_poly_assign(&params.ringq, 0.0, self.e_msg_stddev1, p)
        });
        let mut e2 = vec![params.ringq.new_ntt_poly(); params.l];
        e2.iter_mut().for_each(|p| {
            self.sampler
                .sample_poly_assign(&params.ringq, 0.0, self.e_msg_stddev2, p)
        });

        // Line 13
        let mut eps0 = vec![params.ringq.new_ntt_poly(); params.munu];
        eps0.iter_mut().for_each(|p| {
            self.sampler
                .sample_poly_assign(&params.ringq, 0.0, self.e_rand_stddev1, p)
        });
        let mut eps1 = vec![params.ringq.new_ntt_poly(); params.munu];
        eps1.iter_mut().for_each(|p| {
            self.sampler
                .sample_poly_assign(&params.ringq, 0.0, self.e_rand_stddev1, p)
        });
        let mut eps2 = vec![params.ringq.new_ntt_poly(); params.munu];
        eps2.iter_mut().for_each(|p| {
            self.sampler
                .sample_poly_assign(&params.ringq, 0.0, self.e_rand_stddev2, p)
        });

        // Line 14
        let ecommit0 = self.comitter.commit(&e0, &eps0);
        let ecommit1 = self.comitter.commit(&e1, &eps1);
        let ecommit2 = self.comitter.commit(&e2, &eps2);

        // Line 15
        ecommit0
            .iter()
            .chain(&ecommit1)
            .chain(&ecommit2)
            .for_each(|p| self.oracle.write_poly(p));

        // Line 29
        // We early sample x
        self.oracle.finalize();
        let x = self.oracle.read_range(params.p);
        self.oracle.write_u256(x);

        // Line 20
        let mut edec0 = vec![U256::from(0); params.l * params.m];
        self.encoder.decode_chunk_assign(&e0, &mut edec0);
        let mut edec1 = vec![U256::from(0); params.l * params.m];
        self.encoder.decode_chunk_assign(&e1, &mut edec1);
        let mut edec2 = vec![U256::from(0); params.l * params.m];
        self.encoder.decode_chunk_assign(&e2, &mut edec2);

        // Line 21
        let mut f = vec![U256::from(0); params.M];
        // B.T*V*A*w
        A.mul_vec_assign(&w, &mut buff_vec);
        for i in 0..params.M {
            buff_vec[i] = (buff_vec[i] * v[i]) % params.p; // V*A*w
        }
        B.transpose_mul_vec_assign(&buff_vec, &mut f);

        // Line 22
        let mut g0 = U256::from(0);
        for i in 0..params.M {
            let mut avb = (v[i] * b[i]) % params.p;
            avb = (avb * a[i]) % params.p;
            let cv = (c[i] * v[i]) % params.p;

            if avb >= cv {
                g0 += avb - cv;
            } else {
                g0 += params.p - cv + avb;
            }
        }
        g0 %= params.p;
        let mut g0poly = params.ringmul.new_poly();
        params
            .ringmul
            .set_coeff(&mut g0poly, params.M + 2 * params.l * params.m, g0);
        params.ringmul.ntt(&mut g0poly);

        let mut g1 = U256::from(0);
        for i in 0..params.M {
            g1 += (d_reuse[i] * w[i]) % params.p;
        }
        g1 %= params.p;
        let mut g1xpoly = params.ringmul.new_poly();
        params.ringmul.set_coeff(
            &mut g1xpoly,
            params.M + 2 * params.l * params.m,
            (x * g1) % params.p,
        );
        params.ringmul.ntt(&mut g1xpoly);

        // Line 23
        let mut tpoly = params.ringmul.new_poly();
        for i in 0..params.M {
            params
                .ringmul
                .set_coeff(&mut tpoly, params.l * params.m + i, tdec[i]);
        }
        params.ringmul.ntt(&mut tpoly);

        let mut dpoly = params.ringmul.new_poly();
        for i in 0..params.M {
            params
                .ringmul
                .set_coeff(&mut dpoly, params.l * params.m + params.M - i, d[i]);
        }
        params.ringmul.ntt(&mut dpoly);

        // Line 24 - 25
        let mut epoly0 = params.ringmul.new_poly();
        for i in 0..params.l * params.m {
            params.ringmul.set_coeff(&mut epoly0, i, edec0[i]);
        }
        params.ringmul.ntt(&mut epoly0);

        let mut epoly1 = params.ringmul.new_poly();
        for i in 0..params.l * params.m {
            params
                .ringmul
                .set_coeff(&mut epoly1, params.l * params.m - 1 - i, edec1[i]);
        }
        params.ringmul.ntt(&mut epoly1);

        let mut epoly2 = params.ringmul.new_poly();
        for i in 0..params.l * params.m {
            params.ringmul.set_coeff(&mut epoly2, i, edec2[i]);
        }
        params.ringmul.ntt(&mut epoly2);

        // Line 26
        let mut fxpoly = params.ringmul.new_poly();
        for i in 0..params.M {
            params.ringmul.set_coeff(
                &mut fxpoly,
                params.l * params.m + params.M - i,
                (x * f[i]) % params.p,
            );
        }
        params.ringmul.ntt(&mut fxpoly);

        let mut wxpolyneg = params.ringmul.new_poly();
        for i in 0..params.M {
            params.ringmul.set_coeff(
                &mut wxpolyneg,
                params.l * params.m + i,
                params.p - ((x * w[i]) % params.p),
            );
        }
        params.ringmul.ntt(&mut wxpolyneg);

        // Line 27
        let e0t = params.ringmul.add(&epoly0, &tpoly);
        let e1d = params.ringmul.add(&epoly1, &dpoly);

        let mut h = params.ringmul.sub(&g0poly, &epoly2);
        params.ringmul.mul_add_inplace(&e0t, &e1d, &mut h);

        let mut h1 = g1xpoly.clone();
        params.ringmul.mul_add_inplace(&e0t, &fxpoly, &mut h1);
        params.ringmul.mul_add_inplace(&e1d, &wxpolyneg, &mut h1);
        // params.ringmul.scalar_mul_inplace(x, &mut h1);
        params.ringmul.add_inplace(&h1, &mut h);

        params.ringmul.intt(&mut h);

        // Line 30
        let mut hdec = vec![U256::from(0); 2 * (params.M + params.l * params.m) + 1]; // Add one auxillary zero
        let pbig = Integer::from(params.p.as_u128());
        for i in 0..2 * (params.M + params.l * params.m) + 1 {
            let tmp = params.ringmul.get_coeff_balanced(&h, i).rem_euc(&pbig);
            hdec[i] = U256::from(tmp.to_u128().unwrap());
        }

        // Line 32 - 33
        let mut henc = vec![vec![params.ringq.new_ntt_poly(); params.l]; 2 * params.k + 2];
        let mut eta = vec![vec![params.ringq.new_ntt_poly(); params.munu]; 2 * params.k + 2];
        let mut hcommit = vec![vec![params.ringq.new_ntt_poly(); params.mu]; 2 * params.k + 2];
        for i in 0..2 * params.k + 2 {
            if i < params.k + 2 {
                self.encoder.encode_randomized_chunk_assign(
                    &hdec[i * params.l * params.m..(i + 1) * params.l * params.m],
                    params.s1,
                    &mut henc[i],
                );
            } else {
                self.encoder.encode_randomized_chunk_assign(
                    &hdec[i * params.l * params.m + 1..(i + 1) * params.l * params.m + 1],
                    params.s1,
                    &mut henc[i],
                );
            }

            eta[i].iter_mut().for_each(|p| {
                self.sampler
                    .sample_poly_assign(&params.ringq, 0.0, params.sig2, p)
            });

            self.comitter
                .commit_assign(&henc[i], &eta[i], &mut hcommit[i]);

            hcommit[i].iter().for_each(|p| self.oracle.write_poly(p));
        }

        // Line 36
        let mut m_ajtai = vec![vec![params.ringq.new_poly(); params.l]; 4 * params.k + 2];
        m_ajtai[0..params.k].clone_from_slice(&t);
        m_ajtai[params.k..2 * params.k].clone_from_slice(&denc);
        m_ajtai[2 * params.k..4 * params.k + 2].clone_from_slice(&henc);
        // Line 37
        let mut mu_ajtai = vec![vec![params.ringq.new_poly(); params.munu]; 4 * params.k + 2];
        mu_ajtai[0..params.k].clone_from_slice(&tau);
        mu_ajtai[params.k..2 * params.k].clone_from_slice(&delta);
        mu_ajtai[2 * params.k..4 * params.k + 2].clone_from_slice(&eta);
        // Line 38
        let mut mcommit_ajtai = vec![vec![params.ringq.new_poly(); params.mu]; 4 * params.k + 2];
        mcommit_ajtai[0..params.k].clone_from_slice(&tcommit);
        mcommit_ajtai[params.k..2 * params.k].clone_from_slice(&dcommit);
        mcommit_ajtai[2 * params.k..4 * params.k + 2].clone_from_slice(&hcommit);

        // Line 39
        let ajtai_proof = self.ajprover.prove(&m_ajtai, &mu_ajtai, &mcommit_ajtai);
        if !self.ajverifier.verify(&mcommit_ajtai, &ajtai_proof) {
            panic!("Ajtai proof verification failed");
        }

        // Line 41
        self.oracle.finalize();
        let y = self.oracle.read_range(params.p);
        self.oracle.write_u256(y);

        // Line 42
        let ylm = mod_exp(y, params.l * params.m, params.p);
        let mut yexp = vec![U256::from(1); 2 * params.k + 2];
        let mut yenc = vec![params.ringq.new_poly(); 2 * params.k + 2];
        // y^(i+1)*lm
        yexp[0] = ylm;
        self.encoder.encode_assign(&[yexp[0]], &mut yenc[0]);
        for i in 1..params.k {
            yexp[i] = (yexp[i - 1] * ylm) % params.p;
            self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
        }

        let mut tout = vec![params.ringq.new_poly(); params.l];
        for i in 0..params.l {
            tout[i].clone_from(&e0[i]);
            for j in 0..params.k {
                params
                    .ringq
                    .mul_add_inplace(&t[j][i], &yenc[j], &mut tout[i]);
            }
        }

        // Line 43
        let mut tauout = vec![params.ringq.new_poly(); params.munu];
        for i in 0..params.munu {
            tauout[i].clone_from(&eps0[i]);
            for j in 0..params.k {
                params
                    .ringq
                    .mul_add_inplace(&tau[j][i], &yenc[j], &mut tauout[i]);
            }
        }

        // Line 44
        // y^M-ilm+1
        yexp[params.k - 1] = (ylm * y) % params.p;
        self.encoder
            .encode_assign(&[yexp[params.k - 1]], &mut yenc[params.k - 1]);
        for i in (0..params.k - 1).rev() {
            yexp[i] = (yexp[i + 1] * ylm) % params.p;
            self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
        }

        let mut dout = vec![params.ringq.new_poly(); params.l];
        for i in 0..params.l {
            dout[i].clone_from(&e1[i]);
            for j in 0..params.k {
                params
                    .ringq
                    .mul_add_inplace(&denc[j][i], &yenc[j], &mut dout[i]);
            }
        }

        // Line 45
        let mut deltaout = vec![params.ringq.new_poly(); params.munu];
        for i in 0..params.munu {
            deltaout[i].clone_from(&eps1[i]);
            for j in 0..params.k {
                params
                    .ringq
                    .mul_add_inplace(&delta[j][i], &yenc[j], &mut deltaout[i]);
            }
        }

        // Line 46
        // y^ilm
        yexp[0] = U256::from(1);
        self.encoder.encode_assign(&[yexp[0]], &mut yenc[0]);
        for i in 1..2 * params.k + 2 {
            yexp[i] = (yexp[i - 1] * ylm) % params.p;
            if i < params.k + 2 {
                self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
            }
        }
        for i in params.k + 2..2 * params.k + 2 {
            yexp[i] = (yexp[i] * y) % params.p;
            self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
        }

        let mut hout = vec![params.ringq.new_poly(); params.l];
        for i in 0..params.l {
            hout[i].clone_from(&e2[i]);
            for j in 0..2 * params.k + 2 {
                params
                    .ringq
                    .mul_add_inplace(&henc[j][i], &yenc[j], &mut hout[i]);
            }
        }

        // Line 47
        let mut etaout = vec![params.ringq.new_poly(); params.munu];
        for i in 0..params.munu {
            etaout[i].clone_from(&eps2[i]);
            for j in 0..2 * params.k + 2 {
                params
                    .ringq
                    .mul_add_inplace(&eta[j][i], &yenc[j], &mut etaout[i]);
            }
        }

        R1CSProof {
            dcommit: dcommit,
            ecommit0: ecommit0,
            ecommit1: ecommit1,
            ecommit2: ecommit2,
            hcommit: hcommit,
            t: tout,
            tau: tauout,
            d: dout,
            delta: deltaout,
            h: hout,
            eta: etaout,
        }
    }
}

pub struct R1CSVerifier<'a> {
    pub params: &'a Parameters,
    pub encoder: Encoder<'a>,
    pub oracle: Oracle,

    pub key: &'a CommitKey,

    pub comitter: Comitter<'a>,

    pub tdbound2: f64,
    pub tdboundinf: f64,

    pub taudeltabound2: f64,
    pub taudeltaboundinf: f64,

    pub hbound2: f64,
    pub hboundinf: f64,

    pub etabound2: f64,
    pub etaboundinf: f64,
}

impl<'a> R1CSVerifier<'a> {
    pub fn new(params: &'a Parameters, key: &'a CommitKey) -> R1CSVerifier<'a> {
        let kf = params.k as f64;
        let kln = kf * (params.l * params.n) as f64;
        let kln2 = (2.0 * kf + 2.0) * (params.l * params.n) as f64;
        let kmunun = kf * (params.munu * params.n) as f64;
        let kmunun2 = (2.0 * kf + 2.0) * (params.munu * params.n) as f64;
        let bf = params.b as f64;
        let kapf = params.kap as f64;

        let tdbound2 = ((bf + 2.0) * kapf) / 2.0
            * ((bf + 1.0) * kf.sqrt() * params.s1 + params.s2)
            * kln.sqrt();
        let tdboundinf = 5.0 * ((bf + 2.0) * kapf) / 2.0
            * ((bf + 1.0) * kf.sqrt() * params.s1 + params.s2)
            * kf.sqrt();

        let taudeltabound2 =
            ((bf + 2.0) * kapf) / 2.0 * (kf.sqrt() * params.sig1 + params.sig2) * kmunun.sqrt();
        let taudeltaboundinf =
            5.0 * ((bf + 2.0) * kapf) / 2.0 * (kf.sqrt() * params.sig1 + params.sig2) * kf.sqrt();

        let hbound2 = ((bf + 2.0) * kapf) / 2.0
            * ((2.0 * kf + 2.0).sqrt() * params.s1 + params.s2)
            * kln2.sqrt();
        let hboundinf = 5.0 * ((bf + 2.0) * kapf) / 2.0
            * ((bf + 1.0) * (2.0 * kf + 2.0).sqrt() * params.s1 + params.s2)
            * (2.0 * kf + 2.0).sqrt();

        let etabound2 = ((bf + 2.0) * kapf) / 2.0
            * ((2.0 * kf + 2.0).sqrt() * params.sig1 + params.sig2)
            * kmunun2.sqrt();
        let etaboundinf = 5.0 * ((bf + 2.0) * kapf) / 2.0
            * ((2.0 * kf + 2.0).sqrt() * params.sig1 + params.sig2)
            * (2.0 * kf + 2.0).sqrt();

        R1CSVerifier {
            params: params,
            encoder: Encoder::new(params),
            oracle: Oracle::new(),

            key: key,

            comitter: Comitter::new(params, key),

            tdbound2: tdbound2,
            tdboundinf: tdboundinf,

            taudeltabound2: taudeltabound2,
            taudeltaboundinf: taudeltaboundinf,

            hbound2: hbound2,
            hboundinf: hboundinf,

            etabound2: etabound2,
            etaboundinf: etaboundinf,
        }
    }

    pub fn verify(
        &mut self,
        A: &SparseMatrix,
        B: &SparseMatrix,
        C: &SparseMatrix,
        a: &[U256],
        b: &[U256],
        c: &[U256],
        tcommit: &[Vec<Poly>],
        proof: &R1CSProof,
    ) -> bool {
        let params = self.params;

        // t
        let (l2, linf) = proof.t.iter().fold((0.0, 0.0), |acc, p| {
            let (l2, linf) = params.ringq.norm(p);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < self.tdbound2 && linf < self.tdboundinf) {
            return false;
        }
        // d
        let (l2, linf) = proof.d.iter().fold((0.0, 0.0), |acc, p| {
            let (l2, linf) = params.ringq.norm(p);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < self.tdbound2 && linf < self.tdboundinf) {
            return false;
        }
        // tau
        let (l2, linf) = proof.tau.iter().fold((0.0, 0.0), |acc, p| {
            let (l2, linf) = params.ringq.norm(p);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < self.taudeltabound2 && linf < self.taudeltaboundinf) {
            return false;
        }
        // delta
        let (l2, linf) = proof.delta.iter().fold((0.0, 0.0), |acc, p| {
            let (l2, linf) = params.ringq.norm(p);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < self.taudeltabound2 && linf < self.taudeltaboundinf) {
            return false;
        }
        // h
        let (l2, linf) = proof.h.iter().fold((0.0, 0.0), |acc, p| {
            let (l2, linf) = params.ringq.norm(p);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < self.hbound2 && linf < self.hboundinf) {
            return false;
        }
        // eta
        let (l2, linf) = proof.eta.iter().fold((0.0, 0.0), |acc, p| {
            let (l2, linf) = params.ringq.norm(p);
            (acc.0 + l2, f64::max(acc.1, linf))
        });
        if !(l2.sqrt() < self.etabound2 && linf < self.etaboundinf) {
            return false;
        }

        // Recollect random oracle
        // v <- A, B, C, a, b, c, tcommit
        // TODO: Incorrect for now, must feed more information about A, B, C...
        A.values
            .iter()
            .chain(&B.values)
            .chain(&C.values)
            .for_each(|&x| self.oracle.write_u256(x));
        a.iter()
            .chain(b)
            .chain(c)
            .for_each(|&x| self.oracle.write_u256(x));
        tcommit
            .iter()
            .flatten()
            .for_each(|p| self.oracle.write_poly(p));
        self.oracle.finalize();

        let mut v = vec![U256::from(1); params.M];
        v[1] = self.oracle.read_range(params.p);
        for i in 2..params.M {
            v[i] = (v[i - 1] * v[1]) % params.p;
        }

        // w <- v, dcommit
        self.oracle.write_u256(v[1]);
        proof
            .dcommit
            .iter()
            .flatten()
            .for_each(|p| self.oracle.write_poly(p));
        self.oracle.finalize();

        let mut w = vec![U256::from(1); params.M];
        w[1] = self.oracle.read_range(params.p);
        for i in 2..params.M {
            w[i] = (w[i - 1] * w[1]) % params.p;
        }

        // x <- w, ecommit0, ecommit1, ecommit2
        self.oracle.write_u256(w[1]);
        proof
            .ecommit0
            .iter()
            .chain(&proof.ecommit1)
            .chain(&proof.ecommit2)
            .for_each(|p| self.oracle.write_poly(p));
        self.oracle.finalize();

        let x = self.oracle.read_range(params.p);

        // y <- x, hcommit
        self.oracle.write_u256(x);
        proof
            .hcommit
            .iter()
            .flatten()
            .for_each(|p| self.oracle.write_poly(p));
        self.oracle.finalize();

        let y = self.oracle.read_range(params.p);

        let mut buff_vec = vec![U256::from(0); params.M];
        let mut f = vec![U256::from(0); params.M];
        A.mul_vec_assign(&w, &mut buff_vec);
        for i in 0..params.M {
            buff_vec[i] = (buff_vec[i] * v[i]) % params.p; // V*A*w
        }
        B.transpose_mul_vec_assign(&buff_vec, &mut f);

        let mut g0 = U256::from(0);
        for i in 0..params.M {
            let mut avb = (v[i] * b[i]) % params.p;
            avb = (avb * a[i]) % params.p;
            let cv = (c[i] * v[i]) % params.p;

            if avb >= cv {
                g0 += avb - cv;
            } else {
                g0 += params.p + avb - cv;
            }
        }
        g0 %= params.p;

        let mut g1 = U256::from(0);
        let mut g1tmp = vec![U256::from(0); params.M];
        for i in 0..params.M {
            buff_vec[i] = (v[i] * a[i]) % params.p;
        }
        B.transpose_mul_vec_assign(&buff_vec, &mut g1tmp);
        for i in 0..params.M {
            buff_vec[i] = (v[i] * b[i]) % params.p;
        }
        A.transpose_mul_vec_add_assign(&buff_vec, &mut g1tmp);
        C.transpose_mul_vec_sub_assign(&v, &mut g1tmp);
        for i in 0..params.M {
            g1 += (g1tmp[i] * w[i]) % params.p;
        }
        g1 %= params.p;

        // Check commitments
        let mut commit_lhs = vec![params.ringq.new_poly(); params.mu];
        let mut commit_rhs = vec![params.ringq.new_poly(); params.mu];

        let ylm = mod_exp(y, params.l * params.m, params.p);
        let mut yexp = vec![U256::from(0); 2 * params.k + 2];
        let mut yenc = vec![params.ringq.new_poly(); 2 * params.k + 2];

        // y^(i+1)*lm
        yexp[0] = ylm;
        self.encoder.encode_assign(&[yexp[0]], &mut yenc[0]);
        for i in 1..params.k {
            yexp[i] = (yexp[i - 1] * ylm) % params.p;
            self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
        }

        self.comitter
            .commit_assign(&proof.t, &proof.tau, &mut commit_lhs);
        for i in 0..params.mu {
            commit_rhs[i].clone_from(&proof.ecommit0[i]);
            for j in 0..params.k {
                params
                    .ringq
                    .mul_add_inplace(&tcommit[j][i], &yenc[j], &mut commit_rhs[i]);
            }

            if !commit_lhs[i].equal(&commit_rhs[i]) {
                return false;
            }
        }

        // Line 44
        // y^M-ilm+1
        yexp[params.k - 1] = (ylm * y) % params.p;
        self.encoder
            .encode_assign(&[yexp[params.k - 1]], &mut yenc[params.k - 1]);
        for i in (0..params.k - 1).rev() {
            yexp[i] = (yexp[i + 1] * ylm) % params.p;
            self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
        }

        self.comitter
            .commit_assign(&proof.d, &proof.delta, &mut commit_lhs);
        for i in 0..params.mu {
            commit_rhs[i].clone_from(&proof.ecommit1[i]);
            for j in 0..params.k {
                params
                    .ringq
                    .mul_add_inplace(&proof.dcommit[j][i], &yenc[j], &mut commit_rhs[i]);
            }

            if !commit_lhs[i].equal(&commit_rhs[i]) {
                return false;
            }
        }

        // Line 46
        // y^ilm
        yexp[0] = U256::from(1);
        self.encoder.encode_assign(&[yexp[0]], &mut yenc[0]);
        for i in 1..2 * params.k + 2 {
            yexp[i] = (yexp[i - 1] * ylm) % params.p;
            if i < params.k + 2 {
                self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
            }
        }
        for i in params.k + 2..2 * params.k + 2 {
            yexp[i] = (yexp[i] * y) % params.p;
            self.encoder.encode_assign(&[yexp[i]], &mut yenc[i]);
        }

        self.comitter
            .commit_assign(&proof.h, &proof.eta, &mut commit_lhs);
        for i in 0..params.mu {
            commit_rhs[i].clone_from(&proof.ecommit2[i]);
            for j in 0..2 * params.k + 2 {
                params
                    .ringq
                    .mul_add_inplace(&proof.hcommit[j][i], &yenc[j], &mut commit_rhs[i]);
            }

            if !commit_lhs[i].equal(&commit_rhs[i]) {
                return false;
            }
        }

        let mut tdec = vec![U256::from(0); params.l * params.m];
        let mut ddec = vec![U256::from(0); params.l * params.m];
        let mut hdec = vec![U256::from(0); params.l * params.m];
        self.encoder.decode_chunk_assign(&proof.t, &mut tdec);
        self.encoder.decode_chunk_assign(&proof.d, &mut ddec);
        self.encoder.decode_chunk_assign(&proof.h, &mut hdec);

        let mut yexp = vec![U256::from(1); params.l * params.m];
        for i in 1..params.l * params.m {
            yexp[i] = (yexp[i - 1] * y) % params.p;
        }

        let mut yexplm = vec![U256::from(1); params.M];
        yexplm[0] = mod_exp(y, params.l * params.m, params.p);
        for i in 1..params.M {
            yexplm[i] = (yexplm[i - 1] * y) % params.p;
        }

        let yh = (0..params.l * params.m)
            .fold(U256::from(0), |acc, i| (acc + yexp[i] * hdec[i]) % params.p);

        let yt = (0..params.l * params.m)
            .fold(U256::from(0), |acc, i| (acc + yexp[i] * tdec[i]) % params.p);
        let yd = (0..params.l * params.m).fold(U256::from(0), |acc, i| {
            (acc + yexp[params.l * params.m - i - 1] * ddec[i]) % params.p
        });

        let fy = (0..params.M).fold(U256::from(0), |acc, i| {
            let tmp = (yexplm[params.M - i - 1] * f[i]) % params.p;
            (acc + tmp * y) % params.p
        });

        let wy = (0..params.M).fold(U256::from(0), |acc, i| (acc + yexplm[i] * w[i]) % params.p);
        let wyneg = params.p - wy;

        let yM2lm = mod_exp(y, params.M + 2 * params.l * params.m, params.p);
        let g0y = (g0 * yM2lm) % params.p;
        let g1y = (g1 * yM2lm) % params.p;

        let t0 = (((yt * yd) % params.p) + g0y) % params.p;
        let ytfy = (yt * fy) % params.p;
        let ydwyneg = (yd * wyneg) % params.p;
        let mut t1 = (ytfy + ydwyneg + g1y) % params.p;
        t1 = (x * t1) % params.p;

        let rhs = (t0 + t1) % params.p;

        return yh == rhs;
    }
}
