use num_complex::{ComplexFloat};

pub fn eq_within_epsilon<F: ComplexFloat>(a: F, b: F, eps: F::Real) -> bool {
    F::abs(a - b) < eps
}
