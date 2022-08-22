//! Second order sections.
use std::cmp::max;
use std::collections::{HashMap};
use num_complex::{Complex};
use num_traits::{Zero};

use crate::cplxreal::cplxreal;
use crate::errors::Error;
use crate::errors::Error::{IllegalArgument, InternalError};
use crate::filter_design::ZPKCoeffs;
use crate::sos::SosPairing::{Minimal, Nearest};

use crate::util::{is_real, real_array_to_cplx};
use crate::zpk2tf::zpk2tf;

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
/// The method to used to combine pairs of poles and zeros into sections.
pub enum SosPairing {
    /// Minimize peak gain.
    Nearest=1,

    /// Minimize peak gain under the constraint that odd-order systems should retain one section as first order.
    KeepOdd=2,

    /// Similar to [SosPairing::KeepOdd], except no additional poles or zeros are introduced.
    Minimal=3
}

type FnConsumeCplx = fn(&mut HashMap<usize, Complex<f64>>, pred: &dyn Fn(Complex<f64>) -> bool) -> Result<Complex<f64>, Error>;

#[derive(Debug, Clone)]
/// Represents a list of cascaded second order sections.
pub struct Sos {
    /// List of second order sections
    pub(crate) sections: Vec<SosCoeffs>,
}

#[derive(Debug, Clone)]
#[repr(C)]
/// Represents a single second order section.
pub struct SosCoeffs {
    pub(crate) b0: f64,
    pub(crate) b1: f64,
    pub(crate) b2: f64,
    pub(crate) a0: f64,
    pub(crate) a1: f64,
    pub(crate) a2: f64,
}

impl SosCoeffs {
    pub(crate) fn from_array(arr: &[f64;6]) -> Self {
        SosCoeffs {
            b0: arr[0],
            b1: arr[1],
            b2: arr[2],
            a0: arr[3],
            a1: arr[4],
            a2: arr[5],
        }
    }
}



impl Sos {
    /// Create an empty Sos struct.
    pub fn new() -> Self {
        Sos { sections: vec![] }
    }

    /// Return the number of sections.
    pub fn num_sections(&self) -> usize {
        self.sections.len()
    }

    fn from_vec(vec: Vec<[f64;6]>) -> Self {
        let mut sections: Vec<SosCoeffs> = vec![];

        for s_array in vec.iter() {
            sections.push(SosCoeffs::from_array(s_array));
        }
        Sos { sections }
    }
 }

impl Default for Sos {
    fn default() -> Self {
        Self::new()
    }
}

/// Convert zero-pole-gain filter parameters to second-order sections.
///
/// # Arguments
///
/// * `ZPKCoeffs` - Zero Pole Gain representation of the sytem.
/// * `pairing` - Optional [pairing](`SosPairing`) strategy.
///               If pairing is None, pairing is set to 'nearest'.
///
/// # Returns
///
/// The requested filter in [`Sos`] format.
///
/// # Notes

/// The algorithm used to convert ZPK to SOS format is designed to
/// minimize errors due to numerical precision issues.
///
/// The pairing algorithm attempts to minimize the peak gain of each biquadratic
/// section. This is done by pairing poles with the nearest zeros, starting
/// with the poles closest to the unit circle: [SosPairing::Minimal] outputs may not be usable for filtering.
///
/// ## Algorithms:
///
///  The algorithm is a literal re-implementation of the scypi implementation
///  [described here](https://github.com/scipy/scipy/blob/cf11a137486db1a83751eca8e592cca11a93147d/scipy/signal/_filter_design.py#L1336]).
///
pub fn zpk2sos(zpk: &ZPKCoeffs, pairing: Option<SosPairing>) -> Result<Sos, Error> {
    let pairing = match pairing {
        None => { Nearest }
        Some(value) => { value }
    };

    let mut z = zpk.z.clone();
    let mut p = zpk.p.clone();

    if z.is_empty() && p.is_empty() {
        return Ok( Sos::from_vec(vec![[zpk.k, 0.0, 0.0, 1.0, 0.0, 0.0]]) );
    }

    let n_sections:usize;

    match pairing {
        SosPairing::Minimal => {
            if p.len() < z.len() {
                return Err(IllegalArgument("for analog zpk2sos conversion, must have len(p)>=len(z)".to_string()) );
            }

            n_sections = (p.len() + 1) / 2;
        }
        _ => {
            // Ensure we have the same number of poles and zeros, and make copies.
            p.append( &mut vec![Complex::<f64>::zero(); max(z.len() as i32 - p.len() as i32, 0) as usize]);
            z.append( &mut vec![Complex::<f64>::zero(); max(p.len() as i32 - z.len() as i32, 0) as usize]);

            n_sections = (max(p.len(), z.len()) + 1) / 2;

            if p.len() % 2 == 1 && pairing == Nearest {
                p.push(Complex::<f64>::zero());
                z.push(Complex::<f64>::zero());
            }

            assert_eq!(p.len(), z.len());
        }
    }

    // Ensure we have complex conjugate pairs
    // (note that _cplxreal only gives us one element of each complex pair):
     let z = {
         let (mut c,r) = cplxreal(z.to_vec(), None)?;

         c.append(&mut real_array_to_cplx(&r));
         c
     };

     let p = {
         let (mut c,r) = cplxreal(p.to_vec(), None)?;

         c.append(&mut real_array_to_cplx(&r));
         c
     };

    // For digital filters "worst" is the closest to the unit circle
    let consume_worst: FnConsumeCplx  = find_and_remove_closest_to_unit_circle;

    // Convert vecs of poles and zeroes to HasMaps to make it easier to remove items.
    let mut p: HashMap<usize, Complex<f64>> = p.iter().enumerate().map(|(i,c)| (i, *c)).collect::<HashMap<_, _>>();
    let mut z: HashMap<usize, Complex<f64>> = z.iter().enumerate().map(|(i,c)| (i, *c)).collect::<HashMap<_, _>>();

    //Construct the system, reversing order so the "worst" are last
    let mut sos = vec![vec![Complex::<f64>::zero(); 6]; n_sections];

    for si in (0..=n_sections-1).rev() {
        //Select the next "worst" pole, p1.
        let p1 = consume_worst(&mut p, &|_x| {true})?;

        let num_remaining_reals = p.iter().filter(|x| x.1.im == 0.0 ).count();

        // Pair the pole p1 with a zero.

        if is_real(p1) && num_remaining_reals == 0 {
            // Special case 1: Last remaining real pole.

            if pairing != Minimal {
                let z1 = find_and_remove_closest(&mut z, NumberType::Real, &p1)?;

               // sos[si] = _single_zpksos([z1, 0], [p1, 0], 1)
                for (i, e) in single_zpksos(&[z1, Complex::<f64>::zero()], &[p1, Complex::<f64>::zero()], 1.0)?.iter().enumerate() {
                    sos[si][i] = *e;
                }
            } else if !z.is_empty() {
                let z1 = find_and_remove_closest(&mut z, NumberType::Real, &p1)?;

                // sos[si] = _single_zpksos([z1], [p1], 1)
                for (i, e) in single_zpksos(&[z1], &[p1], 1.0)?.iter().enumerate() {
                    sos[si][i] = *e;
                }
            } else {
                // sos[si] = _single_zpksos([], [p1], 1)
                for (i, e) in single_zpksos(&[], &[p1], 1.0)?.iter().enumerate() {
                    sos[si][i] = *e;
                }
            }
        }
        else if p.len() + 1 == z.len() && !is_real(p1) &&
            p.iter().filter(|x| is_real(*x.1)).count() == 1 &&
            z.iter().filter(|x| is_real(*x.1)).count() == 1
        {
            // Special case 2: There's one real pole and one real zero
            // left, and an equal number of poles and zeros to pair up.
            // We *must* pair with a complex zero

            let z1 = find_and_remove_closest(&mut z, NumberType::Complex, &p1)?;

            //sos[si] = _single_zpksos([z1, z1.conj()], [p1, p1.conj()], 1)
            for (i, e) in single_zpksos(&[z1, z1.conj()], &[p1, p1.conj()], 1.0)?.iter().enumerate() {
                sos[si][i] = *e;
            }
        }
        else {
            let p2=  match is_real(p1) {
                true => {
                    // Find and consume the worst real pole
                    consume_worst(&mut p, &|x| { is_real(x) } )?
                }
                false => {
                    p1.conj()
                }
            };

            if !z.is_empty() {
                // Find closest zero
                let z1 = find_and_remove_closest(&mut z, NumberType::Any, &p1)?;

                if !is_real(z1) {
                    //sos[si] = _single_zpksos([z1, z1.conj()], [p1, p2], 1)
                    for (i, e) in single_zpksos(&[z1, z1.conj()], &[p1, p2], 1.0)?.iter().enumerate() {
                        sos[si][i] = *e;
                    }
                } else if !z.is_empty() {
                    let z2 = find_and_remove_closest(&mut z, NumberType::Real, &p1)?;

                    //sos[si] = _single_zpksos([z1, z2], [p1, p2], 1)
                    for (i, e) in single_zpksos(&[z1, z2], &[p1, p2], 1.0)?.iter().enumerate() {
                        sos[si][i] = *e;
                    }
                } else {
                    //sos[si] = _single_zpksos([z1], [p1, p2], 1)
                    for (i, e) in single_zpksos(&[z1], &[p1, p2], 1.0)?.iter().enumerate() {
                        sos[si][i] = *e;
                    }
                }

            } else {
                // No more zeros.
                for (i, e) in single_zpksos(&[], &[p1, p2], 1.0)?.iter().enumerate() {
                    sos[si][i] = *e;
                }
            }
        }
    }

    assert!(z.is_empty(), "zpk2sos failed to consume all the zeroes.");
    assert!(p.is_empty(), "zpk2sos failed to consume all the poles.");

    // Put gain in first section.
    for e in sos[0][..3].iter_mut() {
        *e *= zpk.k;
    }

    let mut results:Sos = Sos::new();

    for row in sos.iter() {
        let mut row_array:[f64;6] = [0.0;6];

        for (i,e) in row.iter().enumerate() {
            assert!(e.im == 0.0);
            row_array[i] = e.re;
        }
        results.sections.push( SosCoeffs::from_array(&row_array) );
    }

    return Ok( results );
}

/// Create one second-order section from up to two zeros and pole
fn single_zpksos(z: &[Complex<f64>], p:&[Complex<f64>], k: f64) -> Result<Vec<Complex<f64>>, Error> {
    let mut sos = vec![Complex::<f64>::zero(); 6];
    let ba = zpk2tf(z, p, k)?;

    sos[3-ba.b.len()..3].copy_from_slice(&ba.b);
    sos[6-ba.a.len()..6].copy_from_slice(&ba.a);

    return Ok( sos );
}

fn find_and_remove_closest_to_unit_circle(roots: &mut HashMap<usize, Complex<f64>>,  pred: &dyn Fn(Complex<f64>) -> bool) -> Result< Complex<f64>, Error> {
    if roots.is_empty() {
        return Err(InternalError("roots list should not be empty.".to_string()));
    }

    let mut min_value = f64::MAX;
    let mut key: usize = usize::MAX;

    for (k, root) in roots.iter(){

        if ! pred(*root) {
            continue
        }

        let dist_from_unit_circle = (1.0 - root.norm()).abs();
        if dist_from_unit_circle < min_value {
            key = *k;
            min_value = dist_from_unit_circle;
        }
    }

    match roots.remove(&key) {
        None => { Err(InternalError("Root was not in map.".to_string())) }
        Some(value) => { Ok( value ) }
    }
}

enum NumberType {
    Any,
    Real,
    Complex,
}

fn find_and_remove_closest(roots: &mut HashMap<usize, Complex<f64>>, kind: NumberType, to: &Complex<f64> ) -> Result< Complex<f64>, Error> {
    if roots.is_empty() {
        return Err(InternalError("roots list should not be empty.".to_string()));
    }

    let mut min_value = f64::MAX;
    let mut key: usize = usize::MAX;

    for (k, root) in roots.iter() {
        match kind {
            NumberType::Real => { if root.im != 0.0 {continue}  }
            NumberType::Complex => { if root.im == 0.0 {continue} }
            NumberType::Any => {}
        }

        let dist = ((root - to).norm()).abs();
        if dist < min_value {
            key = *k;
            min_value = dist;
        }
    }

    match roots.remove(&key) {
        None => { Err(InternalError("Root was not in map.".to_string())) }
        Some(value) => { Ok( value ) }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use std::f64::consts::PI;
    use num_complex::{Complex};

    use crate::vec_cplx;
    use crate::assert_approx_eq;
    use crate::filter_design::{butter, FilterType, ZPKCoeffs};

    use crate::sos::{single_zpksos, Sos, SosPairing, zpk2sos};
    use crate::sos::SosPairing::{KeepOdd, Nearest, Minimal};

    const EPS: f64 = f64::EPSILON;

    // Cases that match Gnu Octave
    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_1(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![ (-1.0, 0.0), (-1.0, 0.0) ],
            p: vec_cplx![ (0.57149, 0.29360), (0.57149, -0.29360) ],
            k: 1.0,
        };
        
        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos = Sos::from_vec(vec![[1.0, 2.0, 1.0, 1.0, -1.14298, 0.41280]]); // Octave & MATLAB

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_2(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![ (0.0, 1.0), (0.0, -1.0) ],
            p: vec_cplx![ (0.9, 0.0), (-0.9, 0.0), (0.0, 0.7), (0.0, -0.7) ],
            k: 1.0,
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, 0.0, 1.0, 1.0, 0.0, 0.49],
                 [1.0, 0.0, 0.0, 1.0, 0.0, -0.81]]); // Octave

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_3(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![],
            p: vec_cplx![ (0.8, 0.0), (-0.5, 0.25), (-0.5, -0.25) ],
            k: 1.0,
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, 0.0, 0.0, 1.0, 1.0, 0.3125],
                 [1.0, 0.0, 0.0, 1.0, -0.8, 0.0]]); // Octave, Matlab fails

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_4(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(1.0, 0.0), (1.0, 0.0), (0.0, 0.9), (0.0, -0.9) ],
            p: vec_cplx![ (0.99, 0.01), (0.99, -0.01), (0.1, 0.9), (0.1, -0.9) ],
            k: 1.0,
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, 0.0, 0.81, 1.0, -0.2, 0.82],
                 [1.0, -2.0, 1.0, 1.0, -1.98, 0.9802]]); // Octave, Matlab fails

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_5(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(0.9, 0.1), (0.9, -0.1), (-0.9, 0.0) ],
            p: vec_cplx![ (0.75, 0.25), (0.75, -0.25), (0.9, 0.0) ],
            k: 1.0
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");

        match pairing.unwrap() {
            Nearest => {
                let expected_sos =
                    Sos::from_vec(vec![[1.0, 0.9, 0.0, 1.0, -1.5, 0.625],
                                       [1.0, -1.8, 0.82, 1.0, -0.9, 0.0]]); // Octave, Matlab fails

                sos.assert_approx_equal_to(&expected_sos, 1E-4);
            }
            KeepOdd => {
                let expected_sos =
                    Sos::from_vec(vec![[1.0, -1.8, 0.82, 1.0, -1.5, 0.625],
                                       [1.0, 0.9, 0.0, 1.0, -0.9, 0.0]]); // Octave, Matlab fails

                sos.assert_approx_equal_to(&expected_sos, 1E-4);
            }
            _ => { panic!("Unsuported pairing in unit test."); }
        }

    }

    // Cases that differ from Gnu Octave

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_6(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(-0.3090, 0.9511), (-0.3090, -0.9511), (0.8090, 0.5878),
            (0.8090, -0.5878), (-1.0000, 0.0000)],
            p: vec_cplx![(-0.3026, 0.9312), (-0.3026, -0.9312), (0.7922, 0.5755),
            (0.7922, -0.5755), (-0.9791, 0.0000)],
            k: 1.0,
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, 1.0, 0.0, 1.0, 0.97915, 0.0],
                 [1.0, 0.61803, 1.0, 1.0, 0.60515, 0.95873],
                 [1.0, -1.61803, 1.0, 1.0, -1.58430, 0.95873]]); // Octave, Matlab fails

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_7(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(-1.0, - 1.4142), (-1.0, 1.4142),
            (-0.625, -1.0533), (-0.625, 1.0533)],
            p: vec_cplx![(-0.2, -0.6782), (-0.2, 0.6782),
            (-0.1, -0.5385), (-0.1, 0.5385)],
            k: 4.0,
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[4.0, 8.0, 12.0, 1.0, 0.2, 0.3],
                [1.0, 1.25, 1.5, 1.0, 0.4, 0.5]]); //  Matlab

        sos.assert_approx_equal_to(&expected_sos, 1E-3);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_basic_8(#[case] pairing:Option<SosPairing>) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![],
            p: vec_cplx![(0.2, 0.0), (-0.5, 0.25), (-0.5, -0.25)],
            k: 1.0
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, 0.0, 0.0, 1.0, -0.2, 0.0],
                               [1.0, 0.0, 0.0, 1.0, 1.0, 0.3125]]);

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    /*
      The next two examples are adapted from Leland B. Jackson,
      "Digital Filters and Signal Processing (1995) p.400:
      http://books.google.com/books?id=VZ8uabI1pNMC&lpg=PA400&ots=gRD9pi8Jua&dq=Pole%2Fzero%20pairing%20for%20minimum%20roundoff%20noise%20in%20BSF.&pg=PA400#v=onepage&q=Pole%2Fzero%20pairing%20for%20minimum%20roundoff%20noise%20in%20BSF.&f=false
     */

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_jackson_1(#[case] pairing:Option<SosPairing>) {
    let deg2rad = PI / 180.0;
        let k = 1.0;

        let thetas = [22.5, 45.0, 77.5];
        let mags = [0.8, 0.6, 0.9];

        let i = Complex::<f64>::new(0.0, 1.0);

        let mut z =  thetas
                .iter()
                .map(|theta| (theta * deg2rad * i).exp() )
                .collect::<Vec<Complex<f64>>>();
        z.append(&mut z.iter().map(|x| x.conj()).collect::<Vec<Complex<f64>>>());

        let mut p = thetas.iter().zip(mags).map(|(theta, mag)| mag * (theta * deg2rad * i).exp() ).collect::<Vec<Complex<f64>>>();
        p.append(&mut p.iter().map(|x| x.conj()).collect::<Vec<Complex<f64>>>());

        let sos = zpk2sos(&ZPKCoeffs {z, p, k}, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, -1.41421, 1.0, 1.0, -0.84853, 0.36],
                               [1.0, -1.84776, 1.0, 1.0, -1.47821, 0.64],
                               [1.0, -0.43288, 1.0, 1.0, -0.38959, 0.81]]);

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[rstest]
    #[case(Some(Nearest))]
    #[case(Some(KeepOdd))]
    fn test_jackson_2(#[case] pairing:Option<SosPairing>) {
        let deg2rad = PI / 180.0;
        let k = 1.0;

        let thetas_z = [85.0, 10.0];
        let thetas_p = [22.5, 45.0, 77.5];

        let mags = [0.8, 0.6, 0.9];

        let i = Complex::<f64>::new(0.0, 1.0);

        let mut z =  thetas_z
            .iter()
            .map(|theta| (theta * deg2rad * i).exp() )
            .collect::<Vec<Complex<f64>>>();
        z.append(&mut z.iter().map(|x| x.conj()).collect::<Vec<Complex<f64>>>());
        z.append(&mut vec_cplx!((1.0, 0.0), (-1.0, 0.0)));

        let mut p = thetas_p.iter().zip(mags).map(|(theta, mag)| mag * (theta * deg2rad * i).exp() ).collect::<Vec<Complex<f64>>>();
        p.append(&mut p.iter().map(|x| x.conj()).collect::<Vec<Complex<f64>>>());

        let sos = zpk2sos(&ZPKCoeffs {z, p, k}, pairing).expect("Call to zpk2sos failed.");
        let expected_sos =
            Sos::from_vec(vec![[1.0, 0.0, -1.0, 1.0, -0.84853, 0.36],
                               [1.0, -1.96962, 1.0, 1.0, -1.47821, 0.64],
                               [1.0, -0.17431, 1.0, 1.0, -0.38959, 0.81]]);

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    // Test Pairing

    #[rstest]
    #[case(Some(Nearest), Sos::from_vec(vec![[1.0, 1.0, 0.5, 1.0, -0.75, 0.0], [1.0, 1.0, 0.0, 1.0, -1.6, 0.65]]))]
    #[case(Some(KeepOdd), Sos::from_vec(vec![[1.0, 1.0, 0.0, 1.0, -0.75, 0.0], [1.0, 1.0, 0.5, 1.0, -1.6, 0.65]]))]
    #[case(Some(Minimal), Sos::from_vec(vec![[0.0, 1.0, 1.0, 0.0, 1.0, -0.75], [1.0, 1.0, 0.5, 1.0, -1.6, 0.65]]))]
    fn test_pairing(#[case] pairing:Option<SosPairing>, #[case] expected_sos: Sos) {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(-1.0, 0.0), (-0.5, -0.5), (-0.5, 0.5)],
            p: vec_cplx![(0.75, 0.0), (0.8, 0.1), (0.8, -0.1)],
            k: 1.0,
        };

        let sos = zpk2sos(&zpk_input, pairing).expect("Call to zpk2sos failed.");

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[test]
    fn test_sos_butter_low_pass() {
        // This test is on the following on python example:
        //    sos = butter(5, [12], 'low', fs=81.0, output='sos')
        let order = 5;
        let cutoff = 12.0;
        let fs = 81.0;

        let zpk_out = butter(order, FilterType::LowPass(cutoff),fs).expect("butter failed.");
        let sos =zpk2sos(&zpk_out, None).expect("zpk2sos failed");

        let expected_sos = Sos::from_vec(vec![
            [ 0.00659168, 0.01318336, 0.00659168, 1.0, -0.33136391, 0.0],
            [ 1.0, 2.0, 1.0, 1.0, -0.72429772, 0.21290681],
            [ 1.0, 1.0, 0.0, 1.0, -0.95708485, 0.60273144]]);

        sos.assert_approx_equal_to(&expected_sos, 1E-4);
    }

    #[test]
    fn test_single_zpksos() {
        let z = vec_cplx![(-1.0, 0.0), (-1.0, 0.0)];
        let p = vec_cplx![(0.57149, 0.2936), (0.57149, -0.2936)];
        let k = 1.0;

        let sos = single_zpksos(&z, &p, k).expect("Call to single_zpksos failed");

        let expected_sos = vec_cplx![(1.0, 0.0), (2.0, 0.0), (1.0, 0.0), (1.0, 0.0), (-1.14298, 0.0), (0.41280178, 0.0)];

        assert_eq!(sos.len(), expected_sos.len());

        for (i, _) in sos.iter().enumerate() {
            assert_approx_eq!(sos[i], expected_sos[i], 1E-6);
        }
    }
}