use ark_bn254::{Fr, G1Projective};
use ark_std::UniformRand;

use ark_std::Zero;
use pcs::common::polynomial::Polynomial;
use pcs::ipa::commit::batch_commit;
use pcs::ipa::verify::compute_s;
use pcs::ipa::{
    commit::commit,
    open::{batch_evaluation_proof, evaluation_proof},
    setup::GlobalIpaParams,
    verify::{batch_verify, verify},
};
use rand::thread_rng;

mod helpers;
use helpers::*;

#[test]
fn test_ipa_proof() {
    let mut rng = thread_rng();
    let degree = 127;

    let poly = Polynomial::<Fr>::from_random_coefficients(degree);
    let global_params = GlobalIpaParams::<G1Projective>::new(degree);

    let point_x = Fr::rand(&mut rng);

    let poly_commitment = commit(&global_params, &poly).expect("Error commiting Polynomial");

    let (a_m, _g_m, l_r_group, f_x, u_values, u_group_element) =
        evaluation_proof(&global_params, &poly, &point_x)
            .expect("Error evaluatiing polynomial proof");

    let result = verify(
        &global_params,
        &poly_commitment,
        &f_x,
        &point_x,
        &l_r_group,
        &a_m,
        &None,
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

    let mut batched_polys_data = generate_polynomials(n_polys, degree, n_omegas);

    batched_polys_data.commitments = batched_polys_data
        .polys
        .iter()
        .map(|p| commit(&global_params, p).unwrap())
        .collect();

    let (q_commit, q_poly, rho) = batch_commit(
        &global_params,
        &batched_polys_data.polys,
        &batched_polys_data.z_poly,
        &batched_polys_data.zi_polys,
    )
    .unwrap();

    let point_x = Fr::rand(&mut rng);
    let (a_m, _g_m, l_r_group, u_values, u_group_element) = batch_evaluation_proof(
        &global_params,
        &batched_polys_data.polys,
        &q_poly,
        &batched_polys_data.z_poly,
        &batched_polys_data.zi_polys,
        &rho[..],
        &point_x,
    )
    .expect("Error evaluatiing batch polynomial proof");

    let result = batch_verify(
        &global_params,
        &batched_polys_data.commitments,
        &q_commit,
        &batched_polys_data.z_poly,
        &batched_polys_data.zi_polys,
        &point_x,
        &Fr::zero(),
        &l_r_group,
        &a_m,
        &None,
        &u_values,
        &rho[..],
        &u_group_element,
    );
    assert!(result, "Polynomial batch commitment verification failed");
}

#[test]
fn test_ipa_split_ivc() {
    let mut rng = thread_rng();
    let n_polys = 200;
    let n_omegas = 20;
    let degree = 127 - n_omegas;
    let n_iterations = 1;

    let global_params = GlobalIpaParams::<G1Projective>::new(4 * degree);
    let mut proofs: Vec<SplitIvcIpaAccumulatorInput> = Vec::with_capacity(n_iterations);

    let mut batched_polys_data = generate_polynomials(n_polys, degree, n_omegas);

    batched_polys_data.commitments = batched_polys_data
        .polys
        .iter()
        .map(|p| commit(&global_params, p).unwrap())
        .collect();

    let (q_commit, q_poly, rho) = batch_commit(
        &global_params,
        &batched_polys_data.polys,
        &batched_polys_data.z_poly,
        &batched_polys_data.zi_polys,
    )
    .unwrap();

    let x_value = Fr::rand(&mut rng);
    let (a_m, g_m, l_r_group, u_values, u_group) = batch_evaluation_proof(
        &global_params,
        &batched_polys_data.polys,
        &q_poly,
        &batched_polys_data.z_poly,
        &batched_polys_data.zi_polys,
        &rho[..],
        &x_value,
    )
    .expect("Error evaluatiing batch polynomial proof");

    let proof = SplitIvcIpaAccumulatorInput {
        batched_polys_data: Some(batched_polys_data),
        commitment: q_commit,
        x_value,
        f_x: Fr::zero(),
        l_r_group,
        g_m,
        a_m,
        u_group,
        u_values,
        rho_values: rho.clone(),
    };

    for _ in 0..n_iterations {
        proofs.push(proof.clone());
    }

    // Accumulator phase

    let zero = Fr::zero();
    let mut accumulator_inputs: [Option<SplitIvcIpaAccumulatorInput>; 2] = [None, None];
    for i in 0..=n_iterations {
        accumulator_inputs[0] = proofs.get(i).cloned();
        let alpha = Fr::rand(&mut rng);
        let acc_x_value = Fr::rand(&mut rng);
        let mut acc_commitment = G1Projective::zero();
        let mut acc_s: Vec<Fr> = vec![];
        // succint check on input[0] and [1]
        if let Some(input) = &accumulator_inputs[0] {
            let batched_polys_data = input.batched_polys_data.as_ref().unwrap();
            let succint_verification1 = batch_verify(
                &global_params,
                &batched_polys_data.commitments,
                &input.commitment,
                &batched_polys_data.z_poly,
                &batched_polys_data.zi_polys,
                &input.x_value,
                &input.f_x,
                &input.l_r_group,
                &input.a_m,
                &Some(input.g_m),
                &input.u_values,
                &input.rho_values,
                &input.u_group,
            );
            assert!(succint_verification1);
            acc_commitment = input.g_m;
            acc_s = compute_s(&input.u_values, 1 << input.l_r_group.len());
        }
        if let Some(input) = &accumulator_inputs[1] {
            let succint_verification2 = verify(
                &global_params,
                &input.commitment,
                &input.f_x,
                &input.x_value,
                &input.l_r_group,
                &input.a_m,
                &Some(input.g_m),
                &input.u_values,
                &input.u_group,
            );
            assert!(succint_verification2);
            acc_commitment = acc_commitment + input.g_m * alpha;
            let s = compute_s(&input.u_values, 1 << input.l_r_group.len());
            acc_s = s
                .iter()
                .enumerate()
                .map(|(idx, s1)| *s1 * alpha + acc_s.get(idx).unwrap_or_else(|| &zero))
                .collect();
        }
        let poly_s = Polynomial::<Fr>::from_vector_coefficients(acc_s);
        let acc_v_value = poly_s.evaluate(&acc_x_value);

        let (acc_a_m, acc_g_m, acc_l_r_group, _f_x, acc_u_values, acc_u_group_element) =
            evaluation_proof(&global_params, &poly_s, &acc_x_value)
                .expect("Error evaluatiing polynomial proof");
        let new_accumulation = SplitIvcIpaAccumulatorInput {
            batched_polys_data: None,
            commitment: acc_commitment,
            x_value: acc_x_value,
            f_x: acc_v_value,
            l_r_group: acc_l_r_group,
            g_m: acc_g_m,
            a_m: acc_a_m,
            u_group: acc_u_group_element,
            u_values: acc_u_values,
            rho_values: rho.clone(),
        };
        accumulator_inputs[1] = Some(new_accumulation);
    }

    // Full verification of accumulator
    let input = accumulator_inputs[1].clone().unwrap();
    let final_verification = verify(
        &global_params,
        &input.commitment,
        &input.f_x,
        &input.x_value,
        &input.l_r_group,
        &input.a_m,
        &None,
        &input.u_values,
        &input.u_group,
    );
    assert!(final_verification);
}
