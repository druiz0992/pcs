use ark_ec::{pairing::Pairing, CurveGroup};
use ark_std::Zero;

use crate::setup::GlobalParams;
use common::polynomial::Polynomial;

/// To commit to a polynomial f(x) = a_0 + a_1 * x + a_2 * x^2 + .... + a_d * x^d, C_f = Sum{i=0,i=d} a_i * [ s_i * G1], where
/// s_i * G is the result of the setup protocol
pub fn commit<P: Pairing>(
    global_params: &GlobalParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
) -> Result<P::G1Affine, String> {
    if polynomial.degree() > global_params.len() {
        return Err(
            "Error committing to Polynomial. Polynomial degree is higher than the number of powers"
                .to_string(),
        )?;
    }

    let mut commitment = P::G1::zero();
    for (i, coeff) in polynomial.coeffs().iter().enumerate() {
        commitment += global_params.g1_get(i).unwrap().to_owned() * coeff;
    }
    Ok(commitment.into_affine())
}

pub fn commit_g2<P: Pairing>(
    global_params: &GlobalParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
) -> Result<P::G2Affine, String> {
    if polynomial.degree() > global_params.len() {
        return Err(
            "Error committing to Polynomial. Polynomial degree is higher than the number of powers"
                .to_string(),
        )?;
    }

    let mut commitment = P::G2::zero();
    for (i, coeff) in polynomial.coeffs().iter().enumerate() {
        commitment += global_params.g2_get(i).unwrap().to_owned() * coeff;
    }
    Ok(commitment.into_affine())
}
