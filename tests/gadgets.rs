use ark_bn254::{Bn254, Fr};
use ark_std::UniformRand;

use pcs::common::{polynomial::Polynomial, utils};
use pcs::kzg::{commit::commit, open::evaluation_proof, setup::GlobalKzgParams, verify::verify};
use rand::thread_rng;

/// proof that q(X) = f(X)/z(X). where Z is a the vanishing polynomial of omega
#[test]
fn test_kzg_zero_test_on_omega() {
    let mut rng = thread_rng();
    let order = 16;
    let q_degree = 100;
    let omega = utils::compute_roots_of_unity::<Fr>(order).unwrap();
    let vanishing_poly = Polynomial::from_monomial_coefficients(omega);
    let global_params = GlobalKzgParams::<Bn254>::new(2 * q_degree);

    let q_poly = Polynomial::from_random_coefficients(q_degree);
    let polynomial = q_poly.clone() * vanishing_poly.clone();

    let point_u = Fr::rand(&mut rng);

    let commit_f = commit(&global_params, &polynomial).unwrap();
    let eval_qu = q_poly.evaluate(&point_u);
    let eval_zu = vanishing_poly.evaluate(&point_u);

    let (proof, eval_fu) = evaluation_proof(&global_params, &polynomial, &point_u)
        .expect("Error evaluatiing polynomial proof");
    let result = verify(&global_params, &commit_f, &proof, &point_u, &eval_fu);

    assert_eq!(eval_fu, eval_qu * eval_zu);
    assert!(result, "Polynomial commitment verification failed");
}
