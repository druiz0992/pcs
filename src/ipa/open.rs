use super::utils;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::Field;
use ark_ff::UniformRand;
use ark_ff::Zero;
use rand::thread_rng;

use super::setup::GlobalIpaParams;
use crate::common::polynomial::Polynomial;

pub fn evaluation_proof<P: CurveGroup + PrimeGroup>(
    global_params: &GlobalIpaParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
    x_value: &P::ScalarField,
) -> Result<
    (
        P::ScalarField, // a[0]
        P,              // G[0] -> needed for accumulator
        Vec<(P, P)>,    // L,R vectors
        P::ScalarField, // evaluation f(x)
        // Below parameters are not necessary in the
        // final protocol as they will be computed by prover and verifier
        // using fiat shamir transform to make protocol non inteactive
        Vec<P::ScalarField>, // u challenges
        P,                   // U group element
    ),
    String,
> {
    let mut coeffs_b = utils::compute_b::<P>(*x_value, polynomial.degree() + 1);
    let mut coeffs_a = polynomial.coeffs().to_vec();
    let mut g_group_elements = global_params.g_coeffs().to_vec();
    let mut n = coeffs_a.len();
    let m = ark_std::log2(n) as usize;
    let mut l_r_group = Vec::with_capacity(ark_std::log2(n) as usize);
    let u_group = utils::compute_u_group_element::<P>();
    let u_values = compute_u_field_values::<P>(m);
    let f_x = polynomial.evaluate(&x_value);

    for j in 0..m {
        let inner_product_a_g =
            utils::inner_product_group::<P>(&coeffs_a[..n / 2], &g_group_elements[n / 2..], n / 2);
        let inner_product_a_b =
            utils::inner_product_field_element::<P>(&coeffs_a[..n / 2], &coeffs_b[n / 2..], n / 2);
        let group_inner_product_a_b = u_group * inner_product_a_b;
        let l_group = inner_product_a_g + group_inner_product_a_b;

        let inner_product_a_g =
            utils::inner_product_group::<P>(&coeffs_a[n / 2..], &g_group_elements[..n / 2], n / 2);
        let inner_product_a_b =
            utils::inner_product_field_element::<P>(&coeffs_a[n / 2..], &coeffs_b[..n / 2], n / 2);
        let group_inner_product_a_b = u_group * inner_product_a_b;
        let r_group = inner_product_a_g + group_inner_product_a_b;
        l_r_group.push((l_group, r_group));

        let u_inverse = u_values[j].inverse().unwrap();
        let u = u_values[j];
        for i in 0..n / 2 {
            coeffs_a[i] = u * coeffs_a[i] + u_inverse * coeffs_a[n / 2 + i];
            coeffs_b[i] = u_inverse * coeffs_b[i] + u * coeffs_b[n / 2 + i];
            g_group_elements[i] = g_group_elements[i] * u_inverse + g_group_elements[n / 2 + i] * u;
        }
        n /= 2;
    }
    Ok((
        coeffs_a[0],
        g_group_elements[0],
        l_r_group,
        f_x,
        u_values,
        u_group,
    ))
}

pub fn batch_evaluation_proof<P: CurveGroup + PrimeGroup>(
    global_params: &GlobalIpaParams<P>,
    polynomials: &[Polynomial<P::ScalarField>],
    q_poly: &Polynomial<P::ScalarField>,
    z_poly: &Polynomial<P::ScalarField>, // z(x) = Product (X-omega), for all omegas in Big Omega
    z_i_poly: &[Polynomial<P::ScalarField>], // z_i(x) = Product (X - omega) for omega_i in Big Omerga - Big Omega i
    rho: &[P::ScalarField],
    x_value: &P::ScalarField,
) -> Result<
    (
        P::ScalarField, // a[0]
        P,              // G[0]
        Vec<(P, P)>,    // L,R vectors
        // Below parameters are not necessary in the
        // final protocol as they will be computed by prover and verifier
        // using fiat shamir transform to make protocol non inteactive
        Vec<P::ScalarField>, // u challenges
        P,                   // U group element
    ),
    String,
> {
    let z_evaluation = z_poly.evaluate(x_value);
    let zi_evaluations: Vec<P::ScalarField> =
        z_i_poly.iter().map(|p| p.evaluate(x_value)).collect();
    let scaled_zi_evaluations: Vec<P::ScalarField> = zi_evaluations
        .iter()
        .enumerate()
        .map(|(idx, v)| *v * rho[idx])
        .collect();

    let g_poly = compute_g_poly(&polynomials, &q_poly, &z_evaluation, &scaled_zi_evaluations);

    let (a_m, g_m, l_r_group, f_x, u_values, u_group_element) =
        evaluation_proof(&global_params, &g_poly, x_value)?;

    if f_x != P::ScalarField::zero() {
        return Err("g_poly should evaluate to zero at point x".to_string());
    }

    Ok((a_m, g_m, l_r_group, u_values, u_group_element))
}

fn compute_u_field_values<P: CurveGroup>(n: usize) -> Vec<P::ScalarField> {
    let mut rng = thread_rng();
    let mut u = Vec::with_capacity(n);
    for _ in 0..n {
        u.push(P::ScalarField::rand(&mut rng));
    }
    u
}

fn compute_g_poly<F: Field>(
    polynomials: &[Polynomial<F>],
    q_poly: &Polynomial<F>,
    z_evaluation: &F,
    scaled_zi_evaluations: &[F],
) -> Polynomial<F> {
    let scaled_q_poly = Polynomial::<F>::from_vector_coefficients(
        q_poly.coeffs().iter().map(|q| *q * z_evaluation).collect(),
    );
    let mut scaled_f_polys = Polynomial::<F>::from_vector_coefficients(vec![]);
    scaled_f_polys = polynomials
        .iter()
        .enumerate()
        .fold(scaled_f_polys, |acc, (idx, p)| {
            let poly_coeffs: Vec<F> = p
                .coeffs()
                .iter()
                .map(|p| *p * scaled_zi_evaluations[idx])
                .collect();
            let poly = Polynomial::<F>::from_vector_coefficients(poly_coeffs);
            poly + acc
        });
    scaled_f_polys - scaled_q_poly
}
