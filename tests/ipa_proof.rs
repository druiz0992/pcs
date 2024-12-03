use ark_bn254::{Fr, G1Projective};
use ark_std::UniformRand;

use pcs::common::polynomial::Polynomial;
use pcs::ipa::{commit::commit, open::evaluation_proof, setup::GlobalIpaParams, verify::verify};
use rand::thread_rng;

#[test]
fn test_ipa_proof() {
    let mut rng = thread_rng();
    let degree = 127;

    let poly = Polynomial::<Fr>::from_random_coefficients(degree);
    let global_params = GlobalIpaParams::<G1Projective>::new(degree);

    let point_x = Fr::rand(&mut rng);

    let poly_commitment = commit(&global_params, &poly).expect("Error commiting Polynomial");

    let (a_m, l_r_group, f_x, u_values, u_group_element) =
        evaluation_proof(&global_params, &poly, &point_x)
            .expect("Error evaluatiing polynomial proof");

    let result = verify(
        &global_params,
        &poly_commitment,
        &f_x,
        &point_x,
        &l_r_group,
        &a_m,
        &u_values,
        &u_group_element,
    );

    assert!(result, "Polynomial commitment verification failed");
}
