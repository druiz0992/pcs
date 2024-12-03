use ark_bn254::{Bn254, Fr};
use ark_std::UniformRand;

use pcs::common::polynomial::Polynomial;
use pcs::kzg::{
    commit::commit,
    open::{batch_evaluation_proof, evaluation_proof},
    setup::GlobalKzgParams,
    verify::{batch_verify, verify},
};
use rand::thread_rng;

#[test]
fn test_kzg_proof() {
    let mut rng = thread_rng();
    let degree = 100;

    let poly = Polynomial::<Fr>::from_random_coefficients(degree);
    let global_params = GlobalKzgParams::<Bn254>::new(degree);

    let point_u = Fr::rand(&mut rng);

    let poly_commitment = commit(&global_params, &poly).expect("Error commiting Polynomial");

    let (proof, eval_u) = evaluation_proof(&global_params, &poly, &point_u)
        .expect("Error evaluatiing polynomial proof");

    let result = verify(&global_params, &poly_commitment, &proof, &point_u, &eval_u);

    assert!(result, "Polynomial commitment verification failed");
}

#[test]
fn test_kzg_batch_proof() {
    let mut rng = thread_rng();
    let degree = 100;
    let n_commits = 10;
    let mut points_u = Vec::with_capacity(n_commits);

    let poly = Polynomial::<Fr>::from_random_coefficients(degree);
    let global_params = GlobalKzgParams::<Bn254>::new(degree);

    for _ in 0..n_commits {
        points_u.push(Fr::rand(&mut rng));
    }

    let poly_commitment = commit(&global_params, &poly).expect("Error commiting Polynomial");

    let (proof, r_poly) = batch_evaluation_proof(&global_params, &poly, &points_u)
        .expect("Error evaluatiing polynomial proof");

    let result =
        batch_verify(&global_params, &poly_commitment, &proof, &points_u, &r_poly).unwrap();

    assert!(result, "Polynomial commitment batch verification failed");
}
