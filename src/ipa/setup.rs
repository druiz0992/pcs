use ark_ec::AffineRepr;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_std::UniformRand;
use rand::thread_rng;

#[derive(Debug, PartialEq)]
pub struct GlobalIpaParams<P: CurveGroup> {
    pub g: Vec<P>,
    pub h: P,
}

impl<P: CurveGroup> GlobalIpaParams<P> {
    pub fn new(max_degree: usize) -> Self {
        ipa_setup(max_degree)
    }

    pub fn len(&self) -> usize {
        self.g.len()
    }

    pub fn g_iter(&self) -> std::slice::Iter<P> {
        self.g.iter()
    }

    pub fn g_get(&self, n: usize) -> Option<&P> {
        self.g.get(n)
    }

    pub fn g_coeffs(&self) -> &[P] {
        &self.g
    }
    pub fn is_empty(&self) -> bool {
        if self.g.len() == 0 {
            return true;
        }
        false
    }

    pub fn h_get(&self) -> P {
        self.h
    }
}

fn ipa_setup<P: CurveGroup>(max_degree: usize) -> GlobalIpaParams<P> {
    let mut rng = thread_rng();

    let r = P::ScalarField::rand(&mut rng);

    let mut g = Vec::with_capacity(max_degree + 1);
    let g1_generator = P::Affine::generator();
    let h = g1_generator * r;

    for _ in 0..=max_degree {
        let s = P::ScalarField::rand(&mut rng);
        g.push(g1_generator * s);
    }

    GlobalIpaParams { g, h }
}
