
use num_complex::{Complex, ComplexFloat};
use crate::cplx;
use crate::errors::Error;
use crate::errors::Error::IllegalArgument;
use crate::util::Keys::{Im, Re};
use crate::util::{sort_lex_with_map};

pub fn cplxreal(mut z: Vec<Complex<f64>>, tol: Option<f64>) -> Result<(Vec<Complex<f64>>, Vec<f64>), Error> {
    if z.is_empty() {
        return Ok((vec![], vec![]));
    }

    let tol = match tol {
        None => { 100.0 * f64::EPSILON }
        Some(value) => { value }
    };

    //let mut z_abs: Vec<Complex<f64>> = z.iter().map(|x| cplx![x.re, x.im.abs()]).collect();
    sort_lex_with_map(&mut z, &[Re, Im], |c: &Complex<f64>| { cplx!(c.re, (c.im).abs() ) }  )?;

    let zr: Vec<f64> = z
        .iter()
        .flat_map(|x| {
            if (x.im).abs() <= tol * x.abs() { Some(x.re) } else { None }
        }).collect();

    if zr.len() == z.len() {
        return Ok((vec![], zr));
    }
    /*
        # Split positive and negative halves of conjugates
        z = z[~real_indices]
        zp = z[z.imag > 0]
        zn = z[z.imag < 0]
    */
    let z: Vec<Complex<f64>> = z
        .iter()
        .flat_map(|x| {
            if (x.im).abs() <= tol * x.abs() { None } else { Some(*x) }
        }).collect();

    let mut zp: Vec<Complex<f64>> = z
        .iter()
        .flat_map(|x| {
            if x.im > 0.0 { Some(*x) } else { None }
        }).collect();

    let mut zn: Vec<Complex<f64>> = z
        .iter()
        .flat_map(|x| {
            if x.im < 0.0 { Some(*x) } else { None }
        }).collect();

    if zp.len() != zn.len() {
        return Err( IllegalArgument("Array contains complex value with no matching conjugate.".to_string()));
    }

    // Find runs of (approximately) the same real part...
    let mut run_start = 0;
    let mut run_end:Option<usize> = None;

    let mut run_indices:Vec<(usize, usize)> = vec![];

    for i in 0..zp.len()-1 {
        let curr = zp[i].re;
        let next = zp[i+1].re;
        let crit = tol * zp[i].abs();

        if next - curr <= crit {
            run_end = Some(i+1);
        } else {

            match run_end {
                None => {}
                Some(run_end) => { run_indices.push((run_start, run_end)); }
            }

            run_start = i+1;
            run_end = None;
        }
    }

    match run_end {
        None => {}
        Some(run_end) => { run_indices.push((run_start, run_end)); }
    }

    // and sort each run by its imaginary part.
    for r in run_indices {
        sort_lex_with_map(&mut zp[r.0..=r.1], &[Im], |c: &Complex<f64>| { cplx!(c.re, (c.im).abs() ) }  )?;
        sort_lex_with_map(&mut zn[r.0..=r.1], &[Im], |c: &Complex<f64>| { cplx!(c.re, (c.im).abs() ) }  )?;
    }

    let mut zc: Vec<Complex<f64>> = vec![];


    for i in 0..zp.len() {
        let p = zp[i];
        let n = zn[i];
        let n_conj = n.conj();

        // Check that negatives match positives
        if (p - n_conj).abs() > tol * n.abs() {
            return Err( IllegalArgument("Array contains complex value with no matching conjugate.".to_string()));
        }

        // Average out numerical inaccuracy in real vs imag parts of pairs
        zc.push((p + n_conj) / 2.0)
    }

    return Ok((zc, zr));
}

#[cfg(test)]
mod tests {
    use num_complex::{Complex};
    use crate::cplxreal::cplxreal;
    use crate::vec_cplx;
    use crate::assert_approx_eq;
    use crate::errors::Error;

    const EPS: f64 = f64::EPSILON;

    #[test]
    fn test_trivial_input() {
        let (zc,zr) = cplxreal(vec![], None).expect("cplxreal failed");
        assert_eq!(zc.len(), 0);
        assert_eq!(zr.len(), 0);

        let (zc,zr) = cplxreal(vec_cplx![(1.0, 0.0)], None).expect("cplxreal failed");
        assert_eq!(zc.len(), 0);
        assert_eq!(zr.len(), 1);
        assert_eq!(zr[0], 1.0);
    }

    #[test]
    fn output_order_1() {
        let a = vec_cplx![
            (0.0, 1.0), (0.0, -1.0), (EPS, 1.0), (EPS,  -1.0), (-EPS, 1.0), (-EPS, -1.0),
            (1.0, 0.0), (4.0, 0.0), (2.0, 0.0), (3.0, 0.0), (0.0, 0.0), (0.0, 0.0),
            (2.0, 3.0), (2.0, -3.0),
            (1.0-EPS, 1.0), (1.0, 2.0), (1.0, -2.0), (1.0+EPS, -1.0),  // sorts out of order
            (3.0, 1.0), (3.0, 1.0), (3.0, 1.0), (3.0, -1.0), (3.0, -1.0), (3.0, -1.0),
            (2.0, -3.0), (2.0, 3.0)];

        let (zc, zr) = cplxreal(a, None).expect("cplxreal failed");

        let zc_expected = vec_cplx![
            (0.0,1.0), (0.0,1.0), (0.0,1.0),
            (1.0,1.0), (1.0,2.0), (2.0,3.0),
            (2.0,3.0), (3.0,1.0), (3.0,1.0),
            (3.0,1.0)
        ];

        let zr_expected = [0.0, 0.0, 1.0, 2.0, 3.0, 4.0];

        assert_eq!(zc.len(), zc_expected.len());

        for (i, _) in zc.iter().enumerate() {
            assert_approx_eq!(zc[i], zc_expected[i], 1E-12);
        }

        assert_eq!(zr.len(), zr_expected.len());

        for (i, _) in zr.iter().enumerate() {
            assert_approx_eq!(zr[i], zr_expected[i], 1E-12);
        }

    }

    #[test]
    fn test_output_order_2() {

        let a = vec_cplx![
            (1.0-EPS, 1.0), (1.0, 2.0), (1.0, -2.0), (1.0+EPS, -1.0), (1.0+EPS, 3.0), (1.0-2.0*EPS, -3.0),
            (0.0, 1.0), (0.0, -1.0), (2.0, 4.0), (2.0, -4.0), (2.0, 3.0), (2.0, -3.0),(3.0, 7.0), (3.0, -7.0), (4.0-EPS, 1.0),
            (4.0+EPS, -2.0), (4.0, -1.0), (4.0-EPS, 2.0)];

        let (zc, zr) = cplxreal(a, None).expect("cplxreal failed");

        let zc_expected = vec_cplx![
            (0.0, 1.0), (1.0, 1.0), (1.0, 2.0),
            (1.0, 3.0), (2.0, 3.0), (2.0, 4.0),
            (3.0, 7.0), (4.0, 1.0), (4.0, 2.0)
        ];

        assert!(zr.is_empty());

        assert_eq!(zc.len(), zc_expected.len());

        for (i, _) in zc.iter().enumerate() {
            assert_approx_eq!(zc[i], zc_expected[i], 1E-12);
        }
    }

    #[test]
    fn test_unmatched_conjugates() {
        // 1+2i is unmatched
        assert_eq!(
            cplxreal(vec_cplx![(1.0, 3.0), (1.0 , -3.0 ), (1.0, 2.0)], None),
            Err( Error::IllegalArgument("Array contains complex value with no matching conjugate.".to_string())));

        // 1+2j and 1-3j are unmatched
        assert_eq!(
            cplxreal(vec_cplx![(1.0, 3.0), (1.0 , -3.0 ), (1.0, 2.0), (1.0, -3.0)], None),
            Err( Error::IllegalArgument("Array contains complex value with no matching conjugate.".to_string())));

        // 1+3j is unmatched
        assert_eq!(
            cplxreal(vec_cplx![(1.0, 3.0), (1.0 , -3.0 ), (1.0, 3.0)], None),
            Err( Error::IllegalArgument("Array contains complex value with no matching conjugate.".to_string())));

        // No pairs
        assert_eq!(
            cplxreal(vec_cplx![(1.0, 3.0)], None),
            Err( Error::IllegalArgument("Array contains complex value with no matching conjugate.".to_string())));
        assert_eq!(
            cplxreal(vec_cplx![(1.0, -3.0)], None),
            Err( Error::IllegalArgument("Array contains complex value with no matching conjugate.".to_string())));

    }
}