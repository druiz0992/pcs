use ark_ec::pairing::Pairing;
use ark_std::One;
use ark_std::Zero;

use super::commit::commit;
use crate::setup::GlobalParams;
use common::polynomial::Polynomial;

/// we want to show that f(u) = v => u is a root of f(x) - v => (X - u) divides f(x) - v => There exists
/// a polynomial q in Fp such that q(x) (x - u) = f(x) - f(u). The evaluation proof process consits on finding q(x) and
/// commitment Cq = q(s) * G
pub fn evaluation_proof<P: Pairing>(
    global_params: &GlobalParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
    u: &P::ScalarField,
) -> Result<(P::G1Affine, P::ScalarField), String> {
    if polynomial.is_zero() {
        return Err("Polynomial is zero".to_string());
    }
    let f_u = polynomial.evaluate(&u);
    let coeff_0 = polynomial.get_coeff(0).unwrap().clone() - f_u;
    let mut numerator_poly = polynomial.clone();
    numerator_poly.set_constant_coeff(coeff_0);

    let denominator_poly = Polynomial::from_vector_coefficients(vec![
        P::ScalarField::zero() - u,
        P::ScalarField::one(),
    ]);
    let (q_poly, r_poly) = Polynomial::<<P as Pairing>::ScalarField>::div_polynomials(
        numerator_poly,
        &denominator_poly,
    )
    .ok_or("Error in polynomial division")?;
    if !r_poly.is_zero() {
        return Err("Poly not divisible".to_string());
    }

    let proof = commit(global_params, &q_poly)?;
    Ok((proof, f_u))
}

pub fn batch_evaluation_proof<P: Pairing>(
    global_params: &GlobalParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
    u: &[P::ScalarField],
) -> Result<(P::G1Affine, Polynomial<P::ScalarField>), String> {
    if polynomial.is_zero() {
        return Err("Polynomial is zero".to_string());
    }
    let monomials = Polynomial::monomial_vector_from_coefficients(u);
    let roots_poly = Polynomial::from_polys(monomials);
    let (_, r_poly) =
        Polynomial::<<P as Pairing>::ScalarField>::div_polynomials(polynomial.clone(), &roots_poly)
            .ok_or("Error in polynomial division")?;
    let numerator_poly = polynomial.clone() - r_poly.clone();
    let (psy_poly, _) =
        Polynomial::<<P as Pairing>::ScalarField>::div_polynomials(numerator_poly, &roots_poly)
            .ok_or("Error in polynomial division")?;

    let proof = commit(global_params, &psy_poly)?;
    Ok((proof, r_poly))
}
