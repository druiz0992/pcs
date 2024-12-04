use ark_ec::CurveGroup;
use ark_ff::UniformRand;

use super::setup::GlobalIpaParams;
use super::utils;
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

pub fn batch_commit<P: CurveGroup>(
    global_params: &GlobalIpaParams<P>,
    polynomials: &[Polynomial<P::ScalarField>],
    z_poly: &Polynomial<P::ScalarField>, // z(x) = Product (X-omega), for all omegas in Big Omega
    z_i_poly: &[Polynomial<P::ScalarField>], // z_i(x) = Product (X - omega) for omega_i in Big Omerga - Big Omega i
) -> Result<(P, Polynomial<P::ScalarField>, Vec<P::ScalarField>), String> {
    let mut rng = rand::thread_rng();
    if polynomials.iter().any(|p| p.degree() > global_params.len()) {
        return Err(
            "Error batch committing Polynomials. Some polynomial degree is higher than the number of powers"
                .to_string(),
        );
    }

    if polynomials.len() != z_i_poly.len() {
        return Err("Error batch commiting Polynomials. Number of polynomials doesnt match with number of z polynomials".to_string());
    }

    // these coefficients come from verifier in reality
    let rho = utils::compute_b::<P>(P::ScalarField::rand(&mut rng), polynomials.len());

    let mut q_poly = Polynomial::<P::ScalarField>::from_vector_coefficients(vec![]);
    q_poly = polynomials
        .iter()
        .enumerate()
        .fold(q_poly, |acc, (idx, p)| {
            let mut poly = p.clone() * z_i_poly[idx].clone();
            let poly_coeffs: Vec<P::ScalarField> =
                poly.coeffs().iter().map(|p| *p * rho[idx]).collect();
            poly = Polynomial::<P::ScalarField>::from_vector_coefficients(poly_coeffs);
            poly + acc
        });
    let (q_poly, r_poly) = Polynomial::<P::ScalarField>::div_polynomials(q_poly, &z_poly)
        .ok_or("Error in polynomial division".to_string())?;

    if !r_poly.is_zero() {
        return Err("quotient polynomial should be divisible by z poly".to_string())?;
    }

    let mut commitment = P::zero();
    for (i, coeff) in q_poly.coeffs().iter().enumerate() {
        commitment += global_params.g_get(i).unwrap().to_owned() * coeff;
    }

    Ok((commitment, q_poly, rho))
}
