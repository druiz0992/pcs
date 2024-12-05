#![allow(dead_code)]

use super::setup::GlobalIpaParams;
use super::utils;
use crate::common::polynomial::Polynomial;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::Field;
#[allow(unused_imports)]
use ark_ff::One;

pub fn verify<P: CurveGroup + PrimeGroup>(
    global_params: &GlobalIpaParams<P>,
    commitment: &P,
    f_x: &P::ScalarField,
    x_value: &P::ScalarField,
    l_r_group: &[(P, P)],
    a_0: &P::ScalarField,
    g_0: &Option<P>,
    u_values: &[P::ScalarField],
    u_group: &P,
) -> bool {
    let m = l_r_group.len();
    let n = 1 << m;
    let mut c = *commitment + *u_group * f_x;
    let s = compute_s(&u_values, n);

    for i in 0..m {
        let u = u_values[i];
        let u_inverse = u.inverse().unwrap();
        let (l_group, r_group) = l_r_group[i];
        c += l_group * (u * u) + r_group * (u_inverse * u_inverse);
    }

    // We can compute b0 as below as well

    // let mut g = Polynomial::<P::ScalarField>::from_vector_coefficients(vec![P::ScalarField::one()]);
    //for i in 0..m {
    // let binomial = compute_binomial(u_values[m - i - 1], 1 << i);
    //g = g * binomial;
    //}
    //let b_0_field = g.evaluate(x_value);

    let b_coeffs = utils::compute_b::<P>(*x_value, n);
    let b_0_field = utils::inner_product_field_element::<P>(&s, &b_coeffs, n);
    let g_0_group = g_0.unwrap_or_else(|| utils::inner_product_group(&s, &global_params.g, n));

    (g_0_group + *u_group * b_0_field) * a_0 == c
}

pub fn batch_verify<P: CurveGroup + PrimeGroup>(
    global_params: &GlobalIpaParams<P>,
    commitments_f: &[P],
    commitment_q: &P,
    z_poly: &Polynomial<P::ScalarField>, // z(x) = Product (X-omega), for all omegas in Big Omega
    z_i_poly: &[Polynomial<P::ScalarField>], // z_i(x) = Product (X - omega) for omega_i in Big Omerga - Big Omega i
    x_value: &P::ScalarField,
    f_x: &P::ScalarField,
    l_r_group: &[(P, P)],
    a_0: &P::ScalarField,
    g_0: &Option<P>,
    u_values: &[P::ScalarField],
    rho: &[P::ScalarField],
    u_group: &P,
) -> bool {
    let commitment_g =
        preprocess_batch_verify(commitments_f, commitment_q, z_poly, z_i_poly, x_value, rho);

    verify(
        global_params,
        &commitment_g,
        &f_x,
        x_value,
        l_r_group,
        a_0,
        g_0,
        u_values,
        u_group,
    )
}

pub fn preprocess_batch_verify<P: CurveGroup + PrimeGroup>(
    commitments_f: &[P],
    commitment_q: &P,
    z_poly: &Polynomial<P::ScalarField>, // z(x) = Product (X-omega), for all omegas in Big Omega
    z_i_poly: &[Polynomial<P::ScalarField>], // z_i(x) = Product (X - omega) for omega_i in Big Omerga - Big Omega i
    x_value: &P::ScalarField,
    rho: &[P::ScalarField],
) -> P {
    let z_evaluation = z_poly.evaluate(x_value);
    let zi_evaluations: Vec<P::ScalarField> =
        z_i_poly.iter().map(|p| p.evaluate(x_value)).collect();
    let scaled_zi_evaluations: Vec<P::ScalarField> = zi_evaluations
        .iter()
        .enumerate()
        .map(|(idx, v)| *v * rho[idx])
        .collect();

    let mut commitment_linear_combination = P::zero();
    commitment_linear_combination = commitments_f
        .iter()
        .enumerate()
        .fold(commitment_linear_combination, |acc, (idx, c)| {
            *c * scaled_zi_evaluations[idx] + acc
        });
    commitment_linear_combination - *commitment_q * z_evaluation
}

fn compute_binomial<F: Field>(coeff: F, degree: usize) -> Polynomial<F> {
    let mut poly_coeffs = vec![F::zero(); degree + 1];
    poly_coeffs[0] = coeff.inverse().unwrap();
    poly_coeffs[degree] = coeff;
    Polynomial::<F>::from_vector_coefficients(poly_coeffs)
}

pub fn compute_s<F: Field>(coeffs: &[F], n: usize) -> Vec<F> {
    let mut s = Vec::with_capacity(n);
    for i in 0..n {
        let mut s_elem = F::one();
        for j in 0..coeffs.len() {
            let mut c = coeffs[coeffs.len() - j - 1];
            if is_inverse(i, j) {
                c = c.inverse().unwrap();
            }
            s_elem *= c;
        }
        s.push(s_elem);
    }
    s
}

fn is_inverse(row_idx: usize, column_idx: usize) -> bool {
    let count = 1 << column_idx;
    if (row_idx / count) % 2 == 0 {
        return true;
    }
    false
}
