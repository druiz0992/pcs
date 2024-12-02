use ark_bn254::{Fr, G1Projective};
use ark_ec::AffineRepr;
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup};
use ark_ff::Field;
use ark_std::{One, UniformRand};
use rand::thread_rng;

#[derive(Debug, PartialEq)]
pub struct GlobalParams<P: Pairing> {
    pub powers_of_g1: Vec<P::G1Affine>,
    pub powers_of_g2: Vec<P::G2Affine>,
}

impl<P: Pairing> GlobalParams<P> {
    pub fn new(max_degree: usize) -> Self {
        kzg_setup(max_degree)
    }

    pub fn len(&self) -> usize {
        self.powers_of_g1.len()
    }

    pub fn g1_iter(&self) -> std::slice::Iter<P::G1Affine> {
        self.powers_of_g1.iter()
    }
    pub fn g2_iter(&self) -> std::slice::Iter<P::G2Affine> {
        self.powers_of_g2.iter()
    }

    pub fn g1_get(&self, n: usize) -> Option<&P::G1Affine> {
        self.powers_of_g1.get(n)
    }

    pub fn g2_get(&self, n: usize) -> Option<&P::G2Affine> {
        self.powers_of_g2.get(n)
    }

    pub fn is_empty(&self) -> bool {
        if self.powers_of_g1.len() == 0 || self.powers_of_g2.len() == 0 {
            return true;
        }
        false
    }
}

fn kzg_setup<P: Pairing>(max_degree: usize) -> GlobalParams<P> {
    let mut rng = thread_rng();

    let s = P::ScalarField::rand(&mut rng);

    let mut powers_of_g1 = Vec::with_capacity(max_degree + 1);
    let mut powers_of_g2 = Vec::with_capacity(max_degree + 1);
    let mut current_power = P::ScalarField::one();
    let g1_generator = P::G1Affine::generator();
    let g2_generator = P::G2Affine::generator();

    for _ in 0..=max_degree {
        powers_of_g1.push((g1_generator * current_power).into_affine());
        powers_of_g2.push((g2_generator * current_power).into_affine());
        current_power *= s;
    }

    GlobalParams {
        powers_of_g1,
        powers_of_g2,
    }
}

pub fn kzg_lagrange_setup(max_degree: usize) -> Vec<G1Projective> {
    let g1_generator = G1Projective::generator();

    let mut denominators = vec![Fr::one(); max_degree + 1];
    for i in 0..=max_degree {
        let xi = Fr::from((i + 1) as u64); // Compute xi as Fr
        for j in 0..=max_degree {
            if i != j {
                let xj = Fr::from((j + 1) as u64); // Compute xj as Fr
                denominators[i] *= xi - xj; // Accumulate denominator product
            }
        }
    }
    let denominators_inverses: Vec<Fr> =
        denominators.iter().map(|d| d.inverse().unwrap()).collect();

    (0..=max_degree)
        .map(|i| g1_generator * denominators_inverses[i])
        .collect()
}

#[cfg(test)]
mod test {
    use ark_bn254::{Bn254, G2Projective};

    use super::*;

    #[test]
    fn test_kzg_setup_length() {
        let degree = 10;
        let global_params = GlobalParams::<Bn254>::new(degree);

        assert_eq!(global_params.len(), degree + 1);
    }

    #[test]
    fn test_kzg_setup_first_g1_element() {
        let degree = 10;
        let global_params = GlobalParams::<Bn254>::new(degree);
        let g1_generator = G1Projective::generator();

        assert_eq!(global_params.g1_iter().next().unwrap(), &g1_generator);
    }

    #[test]
    fn test_kzg_setup_first_g2_element() {
        let degree = 10;
        let global_params = GlobalParams::<Bn254>::new(degree);
        let g2_generator = G2Projective::generator();

        assert_eq!(global_params.g2_iter().next().unwrap(), &g2_generator);
    }
    #[test]
    fn test_kzg_setup_distinctness() {
        let max_degree = 10;
        let global_params = GlobalParams::<Bn254>::new(max_degree);

        for i in 0..global_params.len() {
            for j in i + 1..global_params.len() {
                assert_ne!(
                    global_params.g1_get(i),
                    global_params.g1_get(j),
                    "The elements at index {} and {} should be distinct",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_kzg_lagrange_setup_length() {
        let degree = 10;
        let lagrange_points = kzg_lagrange_setup(degree);

        assert_eq!(lagrange_points.len(), degree + 1);
    }
}
