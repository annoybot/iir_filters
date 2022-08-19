#[macro_export]
/// Creates a [`Vec`] of [Complex](num_complex::Complex) numbers.
///
/// # Arguments
///
/// The numbers are specified as a list of tuples.

/// # Example
///
/// ```
/// use iir_filters::vec_cplx;
/// use num_complex::Complex;
///
/// let real:Vec<Complex<f64>> =  vec_cplx![(-15.0, 0.0), (-10.0, -10.0), (-10.0, 10.0)];
/// ```
#[doc(hidden)]
macro_rules! vec_cplx {
    ($(($re:expr, $im:expr)),*) => {
        vec![$(Complex::new($re, $im)),*]
    }
}

/// A shorthand for declaring a [Complex](num_complex::Complex) number.
///
/// # Arguments
//
/// * `$re:expr` - An expression for the real part of the complex number.
/// * `$im:expr` - An expression for the imaginary part of the complex number.
///
/// # Example
///
/// ```
/// use iir_filters::cplx;
/// use num_complex::Complex;
///
/// let real:Complex<f64> =  cplx!(2.0, 0.0);
/// let im:Complex<f32> =  cplx!(0.0, 2.0);
/// ```
#[macro_export]
#[doc(hidden)]
macro_rules! cplx {
    ($re:expr, $im:expr) => {
        {
           let z = Complex::new($re, $im);
           z
        }
    }
}