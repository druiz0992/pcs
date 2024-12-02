use ark_ec::{pairing::Pairing, AffineRepr};

use crate::{
    commit::{commit, commit_g2},
    setup::GlobalParams,
};
use common::polynomial::Polynomial;

/// The verifier has trusted setup params [s_i * G], commitment to polynomial f C_f, u and v (such that f(u) = v)
/// and the commitment to polynomial q C_q as the proof.
/// Verifier accepts proof if (s - u)C_q = C_f - v * G <=>  (x-u) * q(x) = f(x) - v  (prover side)
pub fn verify<P: Pairing>(
    global_params: &GlobalParams<P>,
    commitment_f: &P::G1Affine,
    proof: &P::G1Affine,
    u: &P::ScalarField,
    v: &P::ScalarField,
) -> bool {
    let g1_generator = P::G1Affine::generator();
    let g2_generator = P::G2Affine::generator();

    if global_params.is_empty() {
        return false;
    }

    // s * G2 - u * G2
    let g2_s = global_params.g2_get(1).unwrap();
    let lhs_g2 = *g2_s - g2_generator * u;
    let lhs = P::pairing(proof, lhs_g2.into());

    // C_f - v * G
    let rhs_g1 = *commitment_f - g1_generator * v;
    let rhs = P::pairing(rhs_g1.into(), g2_generator);

    lhs == rhs
}

pub fn batch_verify<P: Pairing>(
    global_params: &GlobalParams<P>,
    commitment_f: &P::G1Affine,
    proof: &P::G1Affine,
    u: &[P::ScalarField],
    r_poly: &Polynomial<P::ScalarField>,
) -> Result<bool, String> {
    let g2_generator = P::G2Affine::generator();

    if global_params.is_empty() {
        return Err("Global parameters are empty".to_string());
    }

    if u.len() != r_poly.degree() + 1 {
        return Err("Number of points doesnt equal number of evaluations".to_string());
    }

    let monomials = Polynomial::monomial_vector_from_coefficients(u);
    let accumulator_poly = Polynomial::from_polys(monomials);
    let commitment_a = commit_g2(global_params, &accumulator_poly)?;

    let commitment_r = commit(global_params, &r_poly)?;

    // s * G2 - u * G2
    let lhs = P::pairing(proof, commitment_a);

    // C_f - v * G
    let rhs_g1 = *commitment_f - commitment_r;
    let rhs = P::pairing(rhs_g1.into(), g2_generator);

    Ok(lhs == rhs)
}
