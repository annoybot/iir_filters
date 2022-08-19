// Sorting complex numbers and supporting code


use std::collections::HashSet;
use num_complex::Complex;
use num_traits::Float;
use crate::errors::Error;
use crate::errors::Error::IllegalArgument;
use crate::cplx;
use num_traits::{One, Zero};

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum Keys {
    Re,
    Im,
}

pub(crate) fn real_array_to_cplx<F: Float>(reals: &[F]) -> Vec<Complex<F>> {
    let mut c:Vec<Complex<F>> = vec![];

    reals.iter().for_each( |x| {
        c.push( cplx!(*x, F::zero()));
    });

    return c;
}

pub(crate) fn sort_lex<F: Float>(list: &mut [Complex<F>], keys: &[Keys]) -> Result<(), Error>
{
    sort_lex_with_map(list, keys, |a| { *a })
}

#[allow(clippy::unwrap_used)]
pub(crate) fn sort_lex_with_map<F: Float, MAP>(list: &mut [Complex<F>], keys: &[Keys],
                                       m: MAP) -> Result<(), Error>
    where
        MAP: Fn(&Complex<F>) -> Complex<F>,
{
    if keys.is_empty() && keys.len() > 2 {
        return Err(IllegalArgument("Keys must be an non empty array of at most two elements.".to_string()));
    }

    {
        let keySet: HashSet<Keys> = keys.iter().copied().collect();

        if keys.len() != keySet.len() {
            return Err(IllegalArgument("Keys must not repeat..".to_string()));
        }
    }

    // Validate that all the elements in list are normal.
    // This is to ensure that unwrapping the result of partial_cmp(), as we do below, is safe.
    // TODO: Check if total_cmp stabilises and consider using that.
    // See: https://github.com/rust-lang/rust/issues/72599
    for (i, e) in list.iter().enumerate() {
        if !(e.re.is_finite() && e.im.is_finite()) {
            return Err(IllegalArgument(format!("Cannot sort because element at index {} is not normal.", i, )));
        }
    }

    if keys.len() == 1 {
        match keys[0] {
            Keys::Re => { list.sort_by(|  a, b| m(a).re.partial_cmp(&m(b).re).unwrap()); }
            Keys::Im => { list.sort_by(|a, b| m(a).im.partial_cmp(&m(b).im).unwrap()); }
        }
    } else { //keys.len() == 2
        match keys[0] {
            Keys::Re => {
                list.sort_by(
                    |a, b|
                        m(a).re.partial_cmp(&m(b).re).unwrap()
                            .then(m(a).im.partial_cmp(&m(b).im).unwrap()));
            }
            Keys::Im => {
                list.sort_by(
                    |a, b|
                        m(a).im.partial_cmp(&m(b).im).unwrap()
                            .then(m(a).re.partial_cmp(&m(b).re).unwrap()));
            }
        }
    }

    Ok(())
}

/// Find the coefficients of a polynomial with the given sequence of roots.
///
/// Returns the coefficients of the polynomial whose leading coefficient
/// is one for the given sequence of zeros (multiple roots must be included
/// in the sequence as many times as their multiplicity).
///
/// Example: The input [-1.0, 0.0, 1.0] represents the product of the binomials (x-1)(x-0)(x+1).
/// The resulting output: [1.0, 0.0, -1.0, 0.0] represents the expanded polynomial: x³+0x²-x+0
///
/// See: https://stackoverflow.com/a/33594706
pub(crate) fn poly(roots: &[Complex<f64>]) -> Result<Vec<Complex<f64>>, Error> {
    if roots.is_empty() {
        return Err(IllegalArgument("You must supply at least one root.".to_string()));
    }

    let dim = roots.len();
    let mut p = vec![Complex::<f64>::zero(); dim];

    // Start by setting the polynomial p = 1 then...
    p.push(Complex::<f64>::one());

    // Compute p * (x - r) = p * x - p * r  ∀ r
    // ⚠️ p * x is computed by shifting values left ( the rightmost element is set to zero).

    for r in roots {
       for i in 0..p.len() {
           if i == 0 {
               p[i] = p[i+1];
           }
           else if i == p.len()-1 { // last element = 0 - r*p[i]
               p[i] = - r* p[i];
           }
           else {
               p[i] = p[i+1] - r* p[i];
           }
       }
    }

    return Ok(p);
}

pub(crate) fn is_real(x: Complex<f64> ) -> bool {
    // FIXME: Do we need to use an ϵ here?
    return x.im == 0.0;
}


#[cfg(test)]
mod tests {
    use num_complex::{Complex, Complex64};
    use crate::util::Keys::{Im, Re};
    use crate::util::{poly, sort_lex};
    use crate::vec_cplx;
    use crate::errors::Error::IllegalArgument;
    use crate::assert_approx_eq;

    #[test]
    fn test_sorting_re_im() {
        let mut list: [Complex<f64>; 6] = [
            Complex64::new(2.0, 1.0),
            Complex64::new(1.0, 7.0),
            Complex64::new(2.0, 14.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(9.0, 0.02),
            Complex64::new(8.0, 43.0),
        ];

        let expected: [Complex<f64>; 6] = [
            Complex64::new(1.0, 7.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(2.0, 1.0),
            Complex64::new(2.0, 14.0),
            Complex64::new(8.0, 43.0),
            Complex64::new(9.0, 0.02),
        ];

        sort_lex(&mut list, &[Re, Im]).expect("sort_lex crashed.");

        for (el, expected) in list.iter().zip(expected) {
            assert_eq!(*el, expected);
        }
    }

    #[test]
    fn test_sorting_re_only() {
        let mut list: [Complex<f64>; 6] = [
            Complex64::new(2.0, 1.0),
            Complex64::new(1.0, 7.0),
            Complex64::new(2.0, 14.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(9.0, 0.02),
            Complex64::new(8.0, 43.0),
        ];

        let expected: [Complex<f64>; 6] = [
            Complex64::new(1.0, 7.0),
            Complex64::new(2.0, 1.0),
            Complex64::new(2.0, 14.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(8.0, 43.0),
            Complex64::new(9.0, 0.02),
        ];

        sort_lex(&mut list, &[Re]).expect("sort_lex crashed.");

        for (el, expected) in list.iter().zip(expected) {
            assert_eq!(*el, expected);
        }
    }

    #[test]
    fn test_sorting_im_only() {
        let mut list: [Complex<f64>; 6] = [
            Complex64::new(2.0, 1.0),
            Complex64::new(1.0, 7.0),
            Complex64::new(2.0, 14.0),
            Complex64::new(2.0, 0.0),
            Complex64::new(9.0, 0.02),
            Complex64::new(8.0, 43.0),
        ];

        let expected: [Complex<f64>; 6] = [
            Complex64::new(2.0, 0.0),
            Complex64::new(9.0, 0.02),
            Complex64::new(2.0, 1.0),
            Complex64::new(1.0, 7.0),
            Complex64::new(2.0, 14.0),
            Complex64::new(8.0, 43.0),
        ];

        sort_lex(&mut list, &[Im]).expect("sort_lex crashed.");

        for (i, el) in list.iter().enumerate() {
            assert_eq!(el, &expected[i]);
        }
    }

    #[test]
    fn test_poly_1() {
        let roots =  vec_cplx![(-1.0, 0.0), (0.0, 0.0), (1.0, 0.0)];

        let coeff = poly(&roots).expect("Call to poly failed.");
        let expected_coeff = vec_cplx![(1.0, 0.0), (0.0, 0.0), (-1.0, 0.0), (0.0, 0.0)];

        assert_eq!(coeff.len(), expected_coeff.len());

        for (i, _) in coeff.iter().enumerate() {
            assert_approx_eq!(coeff[i], expected_coeff[i], 1E-12);
        }
    }

    #[test]
    fn test_poly_2() {
        let roots = vec_cplx![(2.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)];

        let coeff = poly(&roots).expect("Call to poly failed.");
        let expected_coeff = vec_cplx![(1.0, 0.0), (-11.0, 0.0), (44.0, 0.0), (-76.0, 0.0), (48.0, 0.0)];

        assert_eq!(coeff.len(), expected_coeff.len());

        for (i, _) in coeff.iter().enumerate() {
            assert_approx_eq!(coeff[i], expected_coeff[i], 1E-12);
        }
    }

    #[test]
    fn test_poly_3() {
        let roots =  vec_cplx![(2.0, 0.0)];

        let coeff = poly(&roots).expect("Call to poly failed.");
        let expected_coeff = vec_cplx![(1.0, 0.0), (-2.0, 0.0)];

        assert_eq!(coeff.len(), expected_coeff.len());

        for (i, _) in coeff.iter().enumerate() {
            assert_approx_eq!(coeff[i], expected_coeff[i], 1E-12);
        }
    }

    #[test]
    fn test_poly_4() {
        let roots =  vec_cplx![];

        let result = poly(&roots);
        let expected_err = Err(IllegalArgument("You must supply at least one root.".to_string()));
        assert_eq!(result, expected_err);
    }
}




