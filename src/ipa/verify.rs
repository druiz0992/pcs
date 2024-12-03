use crate::common::polynomial::Polynomial;

use super::setup::GlobalIpaParams;
use super::utils;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::Field;
use ark_ff::One;

pub fn verify<P: CurveGroup + PrimeGroup>(
    global_params: &GlobalIpaParams<P>,
    commitment: &P,
    f_x: &P::ScalarField,
    x_value: &P::ScalarField,
    l_r_group: &[(P, P)],
    a_0: &P::ScalarField,
    u_values: &[P::ScalarField],
    u_group: &P,
) -> bool {
    let m = l_r_group.len();
    let n = 1 << m;
    let mut c = *commitment + *u_group * f_x;
    let mut g = Polynomial::<P::ScalarField>::from_vector_coefficients(vec![P::ScalarField::one()]);
    let s = compute_s(&u_values, n);

    for i in 0..m {
        let u = u_values[i];
        let u_inverse = u.inverse().unwrap();
        let (l_group, r_group) = l_r_group[i];
        c += l_group * (u * u) + r_group * (u_inverse * u_inverse);
        let binomial = compute_monomial(u_values[m - i - 1], 1 << i);
        g = g * binomial;
    }

    let b_0_field = g.evaluate(x_value);
    // We can compute b0 as below as well
    //let b_coeffs = utils::compute_b::<P>(*x_value, n);
    //let b_0_field = utils::inner_product_field_element::<P>(&s, &b_coeffs, n);
    let g_0_group = utils::inner_product_group(&s, &global_params.g, g.degree() + 1);

    (g_0_group + *u_group * b_0_field) * a_0 == c
}

fn compute_monomial<F: Field>(coeff: F, degree: usize) -> Polynomial<F> {
    let mut poly_coeffs = vec![F::zero(); degree + 1];
    poly_coeffs[0] = coeff.inverse().unwrap();
    poly_coeffs[degree] = coeff;
    Polynomial::<F>::from_vector_coefficients(poly_coeffs)
}

fn compute_s<F: Field>(coeffs: &[F], n: usize) -> Vec<F> {
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
