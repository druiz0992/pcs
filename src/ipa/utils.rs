use ark_ec::CurveGroup;
use ark_ff::One;
use ark_ff::Zero;
use rand::thread_rng;

pub fn compute_b<P: CurveGroup>(x_value: P::ScalarField, n: usize) -> Vec<P::ScalarField> {
    let mut b = Vec::with_capacity(n);
    let mut powers_x = x_value * x_value;
    b.push(P::ScalarField::one());
    b.push(x_value);
    for _ in 2..n {
        b.push(powers_x);
        powers_x *= x_value;
    }
    b
}

pub fn compute_u_group_element<P: CurveGroup>() -> P {
    let mut rng = thread_rng();
    P::rand(&mut rng)
}

pub fn inner_product_group<P: CurveGroup>(
    coeffs_a: &[P::ScalarField],
    coeffs_b: &[P],
    n_elems: usize,
) -> P {
    let mut group_elem = P::zero();
    for i in 0..n_elems {
        group_elem += coeffs_b[i] * coeffs_a[i];
    }
    group_elem
}
pub fn inner_product_field_element<P: CurveGroup>(
    coeffs_a: &[P::ScalarField],
    coeffs_b: &[P::ScalarField],
    n_elems: usize,
) -> P::ScalarField {
    let mut field_elem = P::ScalarField::zero();
    for i in 0..n_elems {
        field_elem += coeffs_a[i] * coeffs_b[i];
    }
    field_elem
}
