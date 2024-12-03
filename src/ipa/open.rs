use super::utils;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::Field;
use ark_ff::UniformRand;
use rand::thread_rng;

use super::setup::GlobalIpaParams;
use crate::common::polynomial::Polynomial;

pub fn evaluation_proof<P: CurveGroup + PrimeGroup>(
    global_params: &GlobalIpaParams<P>,
    polynomial: &Polynomial<P::ScalarField>,
    x_value: &P::ScalarField,
) -> Result<
    (
        P::ScalarField,
        Vec<(P, P)>,
        P::ScalarField,
        Vec<P::ScalarField>,
        P,
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
    Ok((coeffs_a[0], l_r_group, f_x, u_values, u_group))
}

fn compute_u_field_values<P: CurveGroup>(n: usize) -> Vec<P::ScalarField> {
    let mut rng = thread_rng();
    let mut u = Vec::with_capacity(n);
    for _ in 0..n {
        u.push(P::ScalarField::rand(&mut rng));
    }
    u
}
