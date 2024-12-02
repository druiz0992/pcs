use ark_ff::{FftField, Field};
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::{univariate::DensePolynomial, Polynomial as ArkPolynomial};
use rand::thread_rng;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F: Field>(DensePolynomial<F>);

impl<F: Field> Polynomial<F> {
    pub fn from_vector_coefficients(mut coeffs: Vec<F>) -> Self {
        while coeffs.last() == Some(&F::zero()) {
            coeffs.pop();
        }
        let p = DensePolynomial::<F> { coeffs };
        Self(p)
    }

    pub fn monomial_from_coefficient(coeff: F) -> Self {
        let p = DensePolynomial::<F> {
            coeffs: vec![F::from(0u32) - coeff, F::from(1u32)],
        };
        Self(p)
    }
    pub fn monomial_vector_from_coefficients(coeffs: &[F]) -> Vec<Polynomial<F>> {
        let monomials: Vec<Polynomial<_>> = coeffs
            .iter()
            .map(|r| Polynomial::from_vector_coefficients(vec![F::zero() - F::from(*r), F::one()]))
            .collect();
        monomials
    }

    pub fn from_random_coefficients(degree: usize) -> Self {
        let mut rng = thread_rng();
        let mut coeffs = Vec::with_capacity(degree + 1);
        for _ in 0..=degree {
            coeffs.push(F::rand(&mut rng))
        }
        Polynomial::from_vector_coefficients(coeffs)
    }

    pub fn add_coefficients(coeffs1: Vec<F>, coeffs2: Vec<F>) -> Vec<F> {
        let poly1 = Polynomial::from_vector_coefficients(coeffs1);
        let poly2 = Polynomial::from_vector_coefficients(coeffs2);
        let new_poly = poly1 + poly2;
        new_poly.coeffs().to_vec()
    }

    pub fn is_zero(&self) -> bool {
        self.0.coeffs.len() == 0
    }

    pub fn set_constant_coeff(&mut self, coeff: F) {
        self.0.coeffs[0] = coeff;
    }

    pub fn degree(&self) -> usize {
        self.0.degree()
    }

    pub fn coeffs(&self) -> &[F] {
        &self.0.coeffs
    }

    pub fn get_coeff(&self, index: usize) -> Option<&F> {
        self.0.coeffs.get(index)
    }

    pub fn evaluate(&self, point: &F) -> F {
        self.0.evaluate(&point)
    }

    pub fn sub_polynomials(a: &Self, b: &Self) -> Self {
        Self(a.0.clone() - b.0.clone())
    }

    pub fn add_polynomials(a: &Self, b: &Self) -> Self {
        Self(a.0.clone() + b.0.clone())
    }

    pub fn div_polynomials(a: Self, b: &Self) -> Option<(Self, Self)> {
        let p_a = DenseOrSparsePolynomial::from(a.0);
        let p_b = DenseOrSparsePolynomial::from(&b.0);
        if let Some(p_div) = p_a.divide_with_q_and_r(&p_b) {
            Some((Self(p_div.0), Self(p_div.1)))
        } else {
            None
        }
    }
}

impl<F: FftField> Polynomial<F> {
    pub fn from_monomial_coefficients(monomials: Vec<F>) -> Self {
        let mut poly = Polynomial(DensePolynomial::<F> {
            coeffs: vec![F::from(1u64)],
        });
        for coeff in monomials {
            let new_monomial = Polynomial::monomial_from_coefficient(coeff);
            poly = poly * new_monomial;
        }
        poly
    }

    pub fn from_polys(polys: Vec<Polynomial<F>>) -> Self {
        if polys.len() == 0 {
            return Polynomial::from_vector_coefficients(vec![]);
        }
        let mut polys = polys.clone();
        let mut new_poly = polys.remove(0);
        for poly in polys {
            new_poly = new_poly * poly;
        }
        new_poly
    }
    pub fn lagrange_interpolation(points: &[F], evals: &[F]) -> Polynomial<F> {
        let mut result_coeffs = vec![F::zero(); points.len()];
        for (i, (&xi, &yi)) in points.iter().zip(evals).enumerate() {
            let mut numerator_poly = Polynomial::from_vector_coefficients(vec![F::one()]);
            let mut denominator = F::one();

            for (j, (&xj, _)) in points.iter().zip(evals).enumerate() {
                if i != j {
                    // Numerator: (x - xj)
                    let new_monomial = Polynomial::monomial_from_coefficient(xj);
                    numerator_poly = numerator_poly * new_monomial;

                    // Denominator: (xi - xj)
                    denominator *= xi - xj;
                }
            }

            let scaled_poly = numerator_poly
                .coeffs()
                .iter()
                .map(|coeff| *coeff * yi * denominator.inverse().unwrap())
                .collect::<Vec<_>>();

            result_coeffs = Self::add_coefficients(result_coeffs, scaled_poly);
        }
        Polynomial::from_vector_coefficients(result_coeffs)
    }
}

impl<F: Field> Sub for Polynomial<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Polynomial::sub_polynomials(&self, &rhs)
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Polynomial::add_polynomials(&self, &rhs)
    }
}

impl<F: FftField> Mul for Polynomial<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

#[cfg(test)]
mod test {
    use ark_bn254::Fr;

    use super::Polynomial;

    #[test]
    fn test_new_from_vector() {
        let coeffs = vec![Fr::from(1_u64), Fr::from(2_u64)];
        let poly = Polynomial::from_vector_coefficients(coeffs.clone());

        assert_eq!(poly.0.coeffs[0], coeffs[0]);
        assert_eq!(poly.0.coeffs[1], coeffs[1]);
    }

    #[test]
    fn test_poly_degree() {
        let coeffs = vec![Fr::from(1_u64), Fr::from(2_u64)];
        let poly = Polynomial::from_vector_coefficients(coeffs.clone());

        assert_eq!(poly.degree(), 1);
    }

    #[test]
    fn test_new_from_vector_last_coeff_0() {
        let coeffs = vec![Fr::from(1_u64), Fr::from(10_u64), Fr::from(0_u64)];
        let poly = Polynomial::from_vector_coefficients(coeffs.clone());

        assert_eq!(poly.degree(), 1)
    }

    #[test]
    fn test_lagrange_interpolation() {
        let degree = 10;
        let points: Vec<_> = (0..=degree).map(|p| Fr::from(p as u64)).collect();
        let poly = Polynomial::<Fr>::from_random_coefficients(degree);
        let eval_points: Vec<_> = points.iter().map(|p| poly.evaluate(p)).collect();
        let p = Polynomial::<Fr>::lagrange_interpolation(&points, &eval_points);

        assert_eq!(p, poly)
    }
}
