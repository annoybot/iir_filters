use num_complex::Complex;
use crate::errors::Error;
use crate::util::{poly, sort_lex};
use crate::util::Keys::{Im, Re};
use crate::vec_cplx;

pub struct BACoeffs {
    pub(crate) b: Vec<Complex<f64>>,
    pub(crate) a: Vec<Complex<f64>>,
}

pub(crate) fn zpk2tf(z: &[Complex<f64>], p:&[Complex<f64>], k: f64) -> Result<BACoeffs, Error> {
    if z.is_empty() && p.is_empty() {
        return Ok( BACoeffs{ b: vec_cplx![(1.0, 0.0)], a: vec_cplx![(1.0, 0.0)] } );

    }

    let mut b: Vec<Complex<f64>> = poly(z)?.iter().map(|x| k*x ).collect();
    let mut a = poly(p)?;

    // Use real output if possible. If the complex roots are all complex conjugates,
    // return real roots.
    if all_conjugates(&b)? {
        b = b.iter()
            .map(|x| Complex::<f64>::new(x.re, 0.0))
            .collect::<Vec<Complex<f64>>>();
    }

    if all_conjugates(&a)? {
        a = a.iter()
            .map(|x| Complex::<f64>::new(x.re, 0.0))
            .collect::<Vec<Complex<f64>>>();
    }

    return Ok( BACoeffs{ b, a } );
}

//Returns true if all the complex numbers are conjugates, false otherwise.
fn all_conjugates(roots: &[Complex<f64>]) -> Result<bool, Error> {

    let mut pos_roots:Vec<Complex<f64>> = roots
        .iter()
        .copied()
        .filter(|x| x.im > 0.0)
        .collect();

    let mut neg_roots:Vec<Complex<f64>> = roots
        .iter()
        .filter(|x| x.im < 0.0)
        .map(|x| x.conj())
        .collect();

    if pos_roots.len() == neg_roots.len() {
        sort_lex(&mut pos_roots, &[Re, Im])?;
        sort_lex(&mut neg_roots, &[Re, Im])?;

        for i in 0..pos_roots.len() {
            if pos_roots[i] != neg_roots[i] {
                return Ok( false );
            }
        }
    }

    return Ok( true );
}

#[cfg(test)]
mod tests {
    use num_complex::Complex;
    use crate::zpk2tf::{BACoeffs, zpk2tf};
    use crate::vec_cplx;

    #[test]
    fn test_identity() {
        let z:Vec<Complex<f64>> = vec![];
        let p:Vec<Complex<f64>> = vec![];
        let k = 1.0;

        let ba_out = zpk2tf(&z, &p, k).expect("call to zpk2tf failed.");

        let ba_expected = BACoeffs {
            b: vec_cplx![(1.0, 0.0)],
            a: vec_cplx![(1.0, 0.0)],
        };

        ba_out.assert_approx_equal_to(&ba_expected, 1E-12);
    }

    #[test]
    fn test_1() {
        let z:Vec<Complex<f64>> = vec_cplx![(-1.0,0.0), (-1.0,0.0)];
        let p:Vec<Complex<f64>> = vec_cplx![(0.57149, 0.2936), (0.57149, -0.2936)];
        let k = 1.0;

        let ba_out = zpk2tf(&z, &p, k).expect("call to zpk2tf failed.");

        let ba_expected = BACoeffs {
            b: vec_cplx![(1.0, 0.0), (2.0, 0.0), (1.0, 0.0)],
            a: vec_cplx![(1.0, 0.0), (-1.14298, 0.0), (0.41280178, 0.0)],
        };

        ba_out.assert_approx_equal_to(&ba_expected, 1E-6);
    }
}