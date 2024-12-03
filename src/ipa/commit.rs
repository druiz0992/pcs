use ark_ec::CurveGroup;

use super::setup::GlobalIpaParams;
use crate::common::polynomial::Polynomial;

pub fn commit<P: CurveGroup>(
    global_params: &GlobalIpaParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
) -> Result<P, String> {
    if polynomial.degree() > global_params.len() {
        return Err(
            "Error committing to Polynomial. Polynomial degree is higher than the number of powers"
                .to_string(),
        )?;
    }

    let mut commitment = P::zero();
    for (i, coeff) in polynomial.coeffs().iter().enumerate() {
        commitment += global_params.g_get(i).unwrap().to_owned() * coeff;
    }

    Ok(commitment)
}
