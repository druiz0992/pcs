use ark_bn254::{Fr, G1Projective};
use ark_std::One;
use ark_std::UniformRand;

use pcs::common::polynomial::Polynomial;
use pcs::ipa::commit::batch_commit;
use pcs::ipa::{
    commit::commit,
    open::{batch_evaluation_proof, evaluation_proof},
    setup::GlobalIpaParams,
    verify::{batch_verify, verify},
};
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

// emulates a plonk proof with multiple polynomials and multiple commitments
#[test]
fn test_ipa_batch_proof() {
    let mut rng = thread_rng();
    let n_polys = 200;
    let n_omegas = 20;
    let degree = 127 - n_omegas;

    let global_params = GlobalIpaParams::<G1Projective>::new(4 * degree);
    let mut polys = Vec::with_capacity(n_polys);
    let mut omegas: Vec<Vec<Fr>> = vec![vec![]];
    let mut zi_poly: Vec<Polynomial<Fr>> = vec![];
    let mut z_poly = Polynomial::<Fr>::from_vector_coefficients(vec![Fr::one()]);

    println!("Preparing input polynomials...");
    for i in 0..n_polys {
        omegas.push((0..n_omegas).map(|_| Fr::rand(&mut rng)).collect());
        zi_poly.push(Polynomial::<Fr>::from_monomial_coefficients(
            omegas[i].clone(),
        ));
        polys.push(Polynomial::<Fr>::from_random_coefficients(degree) * zi_poly[i].clone());
        z_poly = z_poly * zi_poly[i].clone();
    }

    for i in 0..n_polys {
        let t = Polynomial::<Fr>::div_polynomials(z_poly.clone(), &zi_poly[i]).unwrap();
        assert!(t.1.is_zero());
        zi_poly[i] = t.0;
    }

    println!("Commiting to individual polynomials...");
    let poly_commitments: Vec<G1Projective> = polys
        .iter()
        .map(|p| commit(&global_params, p).unwrap())
        .collect();

    println!("Commit  q poly...");
    let (q_commit, q_poly, rho) = batch_commit(&global_params, &polys, &z_poly, &zi_poly).unwrap();

    println!("Batch evaluation...");
    let point_x = Fr::rand(&mut rng);
    let (a_m, l_r_group, u_values, u_group_element) = batch_evaluation_proof(
        &global_params,
        &polys,
        &q_poly,
        &z_poly,
        &zi_poly,
        &rho,
        &point_x,
    )
    .expect("Error evaluatiing batch polynomial proof");

    println!("Batch verification...");
    let result = batch_verify(
        &global_params,
        &poly_commitments,
        &q_commit,
        &z_poly,
        &zi_poly,
        &point_x,
        &l_r_group,
        &a_m,
        &u_values,
        &rho,
        &u_group_element,
    );
    assert!(result, "Polynomial batch commitment verification failed");
}
