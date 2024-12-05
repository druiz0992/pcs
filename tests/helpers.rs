#![allow(dead_code)]

use ark_bn254::Fr;
use ark_bn254::G1Projective;
use ark_std::One;
use ark_std::UniformRand;

use pcs::common::polynomial::Polynomial;
use rand::thread_rng;

#[derive(Debug, Clone)]
pub(crate) struct SplitIvcIpaAccumulatorInput {
    pub batched_polys_data: Option<BatchedPolynomialData>,
    pub commitment: G1Projective,
    pub x_value: Fr,
    pub f_x: Fr,
    pub l_r_group: Vec<(G1Projective, G1Projective)>,
    pub g_m: G1Projective,
    pub a_m: Fr,
    pub u_group: G1Projective,
    pub u_values: Vec<Fr>,
    pub rho_values: Vec<Fr>,
}

#[derive(Debug, Clone)]
pub(crate) struct BatchedPolynomialData {
    pub polys: Vec<Polynomial<Fr>>,
    pub z_poly: Polynomial<Fr>,
    pub zi_polys: Vec<Polynomial<Fr>>,
    pub commitments: Vec<G1Projective>,
}

pub(crate) fn generate_polynomials(
    n_polys: usize,
    degree: usize,
    n_omegas: usize,
) -> BatchedPolynomialData {
    let mut rng = thread_rng();

    let mut polys = Vec::with_capacity(n_polys);
    let mut omegas: Vec<Vec<Fr>> = vec![vec![]];
    let mut zi_polys: Vec<Polynomial<Fr>> = vec![];
    let mut z_poly = Polynomial::<Fr>::from_vector_coefficients(vec![Fr::one()]);

    for i in 0..n_polys {
        omegas.push((0..n_omegas).map(|_| Fr::rand(&mut rng)).collect());
        zi_polys.push(Polynomial::<Fr>::from_monomial_coefficients(
            omegas[i].clone(),
        ));
        polys.push(Polynomial::<Fr>::from_random_coefficients(degree) * zi_polys[i].clone());
        z_poly = z_poly * zi_polys[i].clone();
    }

    for i in 0..n_polys {
        let t = Polynomial::<Fr>::div_polynomials(z_poly.clone(), &zi_polys[i]).unwrap();
        assert!(t.1.is_zero());
        zi_polys[i] = t.0;
    }
    BatchedPolynomialData {
        polys,
        z_poly,
        zi_polys,
        commitments: vec![],
    }
}
