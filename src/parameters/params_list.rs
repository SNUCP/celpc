use super::ParametersLiteral;

impl ParametersLiteral {
    pub fn logn_19_logp_256() -> ParametersLiteral {
        ParametersLiteral {
            ajtai_size: 1,
            ajtai_rand_size: 1 + 2,

            degree: 1 << 19,
            bigint_commit_size: 1 << 12,

            modulus_base: 63388,
            digits: 16,

            ring_degree: 1 << 11,
            ring_modulus: vec![72057594037948417, 72057594037641217],

            commit_std_dev: 10.0,
            opening_proof_std_dev: 34.0,
            blind_std_dev: 5202283.0,

            commit_rand_std_dev: 20.0,
            opening_proof_rand_std_dev: 68.0,
            blind_rand_std_dev: 10404567.0,

            open_proof_bound: f64::exp2(35.7),
            eval_bound: f64::exp2(54.6),
        }
    }

    pub fn logn_21_logp_256() -> ParametersLiteral {
        ParametersLiteral {
            ajtai_size: 1,
            ajtai_rand_size: 1 + 2,

            degree: 1 << 21,
            bigint_commit_size: 1 << 13,

            modulus_base: 63388,
            digits: 16,

            ring_degree: 1 << 11,
            ring_modulus: vec![72057594037948417, 72057594037641217],

            commit_std_dev: 10.0,
            opening_proof_std_dev: 34.0,
            blind_std_dev: 5202283.0,

            commit_rand_std_dev: 20.0,
            opening_proof_rand_std_dev: 68.0,
            blind_rand_std_dev: 10404567.0,

            open_proof_bound: f64::exp2(36.8),
            eval_bound: f64::exp2(55.7),
        }
    }

    pub fn logn_23_logp_256() -> ParametersLiteral {
        ParametersLiteral {
            ajtai_size: 1,
            ajtai_rand_size: 1 + 2,

            degree: 1 << 23,
            bigint_commit_size: 1 << 14,

            modulus_base: 63388,
            digits: 16,

            ring_degree: 1 << 11,
            ring_modulus: vec![72057594037948417, 72057594037641217],

            commit_std_dev: 10.0,
            opening_proof_std_dev: 34.0,
            blind_std_dev: 5202283.0,

            commit_rand_std_dev: 20.0,
            opening_proof_rand_std_dev: 68.0,
            blind_rand_std_dev: 10404567.0,

            open_proof_bound: f64::exp2(38.0),
            eval_bound: f64::exp2(56.9),
        }
    }

    pub fn logn_25_logp_256() -> ParametersLiteral {
        ParametersLiteral {
            ajtai_size: 1,
            ajtai_rand_size: 1 + 2,

            degree: 1 << 25,
            bigint_commit_size: 1 << 15,

            modulus_base: 63388,
            digits: 16,

            ring_degree: 1 << 11,
            ring_modulus: vec![72057594037948417, 72057594037641217],

            commit_std_dev: 10.0,
            opening_proof_std_dev: 34.0,
            blind_std_dev: 5202283.0,

            commit_rand_std_dev: 20.0,
            opening_proof_rand_std_dev: 68.0,
            blind_rand_std_dev: 10404567.0,

            open_proof_bound: f64::exp2(39.2),
            eval_bound: f64::exp2(58.1),
        }
    }
}
