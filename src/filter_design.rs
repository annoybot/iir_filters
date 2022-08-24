//! Functions for designing iir filters.
//!
//! Only Butterworth filters are supoorted at the moment.

use std::f64::consts::PI;
use std::ops::Range;
use num_complex::{Complex, Complex64};
use num_traits::{One, Pow, Zero};
use crate::errors::Error;
use crate::errors::Error::{IllegalArgument, InternalError};
use crate::filter_design::FilterType::{BandPass};
use crate::cplx;


#[derive(Debug, PartialEq)]
/// Represents the type of filter to design, low pass, high pass etc.
///
/// The enum values represent the cutoff frequencies: a single frequency for [LowPass][FilterType::LowPass]  and [HighPass][FilterType::HighPass] ,
/// and two frequencies, high and low respectively for [BandPass][FilterType::BandPass] and [BandStop][FilterType::BandStop].
#[allow(missing_docs)]
pub enum FilterType {
    LowPass(f64),
    HighPass(f64),
    BandPass(f64, f64),
    BandStop(f64, f64),
}

#[derive(Debug, PartialEq)]
/// Represents the algorithm or 'category' of the filter to design, i.e: Butterworth, Chebyshev, etc...
#[allow(missing_docs)]
enum FilterAlgorithm {
    Butterworth
}

/// Filter parameters in Zero, Pole, Gain format.
pub struct ZPKCoeffs {
    ///Zeroes
    pub z: Vec<Complex64>,
    /// Poles
    pub p: Vec<Complex64>,
    /// Gain
    pub k: f64,
}

/// Designs an Nth-order butterworth filter.
///
/// # Arguments
///
/// * `N` - The order of the filter.
/// * `filter_type` - The type for the filter, i.e. low pass, high pass etc...
/// * `fs` - the sampling frequency.
///
/// # Returns
///
/// The requested filter in Zero Pole Gain format.
///
/// # Example
///```rust
/// use iir_filters::filter_design::butter;
/// use iir_filters::filter_design::FilterType;
/// use iir_filters::filter_design::FilterType::BandPass;
///
/// fn main() -> Result<(), Box<dyn std::error::Error>> {
///     let order = 5;
///     let cutoff_low = 1.0;
///     let cutoff_hi= 10.0;
///     let fs = 81.0;
///
///     let zpk = butter(order,FilterType::BandPass(cutoff_low, cutoff_hi), fs)?;
///
///     return Ok( () );
///  }
/// ```
pub fn butter(N: u32, filter_type: FilterType, fs: f64) -> Result<ZPKCoeffs, Error> {
    iirfilter(N, filter_type, FilterAlgorithm::Butterworth, Some( fs ) )
}

fn butter_internal(N: u32, filter_type: FilterType, fs: Option<f64>) -> Result<ZPKCoeffs, Error> {
    iirfilter(N, filter_type, FilterAlgorithm::Butterworth, fs )
}


/// Designs an Nth-order digital filter. Only butterworth filters are supported at the moment.
///
/// # Arguments
///
/// * `N` - The order of the filter.
/// * `filter_type` - The type for the filter, i.e. low pass, high pass etc...
/// * `fs` - the sampling frequency.
///
/// # Returns
///
/// The requested filter in (Zero Pole Gain)[ZPKCoeffs] format.
///
fn iirfilter(N: u32, filter_type:FilterType, algorithm: FilterAlgorithm, fs: Option<f64>) -> Result<ZPKCoeffs, Error> {
    // Convert enum based frequency parameters to python style Wn array for internal use.
    let mut Wn:Vec<f64> = match filter_type {
        FilterType::LowPass(cutoff) | FilterType::HighPass(cutoff) => { vec![cutoff] }
        BandPass(cutoff_lo, cutoff_hi)  | FilterType::BandStop(cutoff_lo, cutoff_hi) => { vec![cutoff_lo, cutoff_hi] }
    };

    // Validate and scale the frequencies.
    for Wn in &mut Wn {
        if *Wn <= 0.0 {
            return Err(IllegalArgument("Frequencies Wn must not be zero.".to_string()))
        }

        if let Some(fs) = fs {
            *Wn = 2.0 * *Wn / fs
        }
    }

    // Validate the scaled frequencies.
    for Wn in &Wn {
        if *Wn <= 0.0 || *Wn >= 1.0 {
            match fs {
                Some(fs) => {
                    return Err( IllegalArgument(format!("Digital filter critical frequencies must be 0 < Wn < fs/2 (fs={} -> fs/2={})", fs, fs/2.0)) );
                }
                None => {
                    return Err( IllegalArgument("Digital filter critical frequencies must be 0 < Wn < 1".to_string()) );
                }
            }
        }
    }

    // ⚠️ Done scaling frequencies. We now deliberately shadow Wn to make it immutable from here on.
    let  Wn= Wn;

    let fs = 2.0;
    let warped:Vec<f64> = Wn.iter().map(|x| 2.0 * fs * f64::tan(PI * x / fs )).collect();

    let zpk =  match algorithm {
        FilterAlgorithm::Butterworth => { butterap(N)? }
    };

    let zpk = match filter_type {
        FilterType::LowPass(_) => {
            assert_eq!(Wn.len(), 1, "Wn.len() != 1. It should contain a single frequency for filter_type {:?}.", filter_type);
            lp2lp_zpk(&zpk,warped[0])?
        }
        FilterType::HighPass(_) => {
            assert_eq!(Wn.len(), 1, "Wn.len() != 1. It should contain a single frequency for this filter_type {:?}.", filter_type);
            lp2hp_zpk(&zpk,warped[0])?
        }
        FilterType::BandPass(_,_)  => {
            assert_eq!(Wn.len(), 2, "Wn.len() != 2. It should contain exactly two frequencies for this filter_type {:?}.", filter_type);
            let bw = warped[1] - warped[0];
            let wo = f64::sqrt(warped[0] * warped[1]);

            lp2bp_zpk(&zpk, wo, bw)?
        }
        FilterType::BandStop(_,_)  => {
            assert_eq!(Wn.len(), 2, "Wn.len() != 2. It should contain exactly two frequencies for this filter_type {:?}.", filter_type);
            let bw = warped[1] - warped[0];
            let wo = f64::sqrt(warped[0] * warped[1]);

            lp2bs_zpk(&zpk, wo, bw)?
        }
    };

    let zpk = bilinear_zpk(&zpk, fs)?;

    Ok( zpk )
}

/// Return [ZPKCoeffs] for analog prototype of Nth-order Butterworth filter.
/// The filter will have an angular (e.g., rad/s) cutoff frequency of 1.
fn butterap(N:u32) -> Result<ZPKCoeffs, Error> {
    let N2:i32 = N as i32;
    let i = Complex::<f64>::i();

    let m:Vec<Complex<f64>> = Range { start: -N2+1, end: N2}
        .step_by(2)
        .map(|x| Complex::<f64>::new(f64::from(x), 0.0))
        .collect();

    //Middle value is 0 to ensure an exactly real pole
    let p:Vec<Complex<f64>> = m.iter().map(|x|  -((i * PI * x / (2.0 * N2 as f64) ).exp())  ).collect();

    let z:Vec<Complex64> = vec![];
    let k = 1.0;

    Ok( ZPKCoeffs { z, p, k} )
}

// Bandpass
fn lp2bp_zpk(zpk: &ZPKCoeffs, wo:f64, bw:f64) -> Result<ZPKCoeffs, Error>
{
    let degree = relative_degree(zpk)?;

    // Scale poles and zeros to desired bandwidth
    let z_lp:Vec<Complex64> = zpk.z.iter().map(|x| x * bw / 2.0).collect();
    let p_lp:Vec<Complex64> = zpk.p.iter().map(|x| x * bw / 2.0).collect();

    let cplx2 =  cplx!( 2.0, 0.0);

    /* Duplicate poles and zeros and shift from baseband to +wo and -wo

       Python code:
       z_bp = concatenate((z_lp + sqrt(z_lp**2 - wo**2),
                    z_lp - sqrt(z_lp**2 - wo**2)))
       p_bp = concatenate((p_lp + sqrt(p_lp**2 - wo**2),
                    p_lp - sqrt(p_lp**2 - wo**2)))
    */
    let mut z_bp:Vec<Complex<f64>> =
        z_lp.iter()
        .map(|x| x + (x.pow(2.0_f64) - wo.pow(2.0_f64) ).sqrt())
            .chain(z_lp.iter()
                .map(|x| x - (x.pow(2.0_f64) - wo.pow(2.0_f64) ).sqrt() ))
            .map(|x| Complex64::one() * x)
            .collect::<Vec<Complex<f64>>>();

    let p_bp:Vec<Complex<f64>> =
        p_lp.iter()
            .map(|x| x + (x.powc(cplx2) - wo.powf(2.0)).sqrt())
            .chain(p_lp.iter()
                .map(|x| x - (x.powc(cplx2) - wo.powf(2.0) ).sqrt() ))
            .collect::<Vec<Complex<f64>>>();

    //Move degree zeros to origin, leaving degree zeros at infinity for BPF
    z_bp.append(&mut vec![Complex::<f64>::zero(); degree]);

    // Cancel out gain change from frequency scaling
    let k_bp = zpk.k * bw.powf(degree as f64);

    // _filter_design.py line 2850
    Ok( ZPKCoeffs { z:z_bp, p: p_bp, k: k_bp} )
}


// Bandstop
fn lp2bs_zpk(zpk: &ZPKCoeffs, wo:f64, bw:f64) -> Result<ZPKCoeffs, Error>
{
    let cplx2 =  cplx!( 2.0, 0.0);
    let degree = relative_degree(zpk)?;

    // Invert to a highpass filter with desired bandwidth
    let z_hp:Vec<Complex64> = zpk.z.iter().map(|x| (bw / 2.0) / x ).collect();
    let p_hp:Vec<Complex64> = zpk.p.iter().map(|x| (bw / 2.0) / x) .collect();

    /* Duplicate poles and zeros and shift from baseband to +wo and -wo

       Python Code
        z_bs = concatenate((z_hp + sqrt(z_hp**2 - wo**2),
                            z_hp - sqrt(z_hp**2 - wo**2)))
        p_bs = concatenate((p_hp + sqrt(p_hp**2 - wo**2),
                            p_hp - sqrt(p_hp**2 - wo**2)))
     */
    let mut z_bs:Vec<Complex<f64>> =
        z_hp.iter()
            .map(|x| x + (x.pow(2.0_f64) - wo.pow(2.0_f64) ).sqrt())
            .chain(z_hp.iter()
                .map(|x| x - (x.pow(2.0_f64) - wo.pow(2.0_f64) ).sqrt() ))
            .map(|x| Complex64::one() * x)
            .collect::<Vec<Complex<f64>>>();

    let p_bs:Vec<Complex<f64>> =
        p_hp.iter()
            .map(|x| x + (x.powc(cplx2) - wo.powf(2.0)).sqrt())
            .chain(p_hp.iter()
                .map(|x| x - (x.powc(cplx2) - wo.powf(2.0) ).sqrt() ))
            .collect::<Vec<Complex<f64>>>();

    //Move any zeros that were at infinity to the center of the stopband
    z_bs.append(&mut vec![Complex::<f64>::i() * wo; degree]);
    z_bs.append(&mut vec![-Complex::<f64>::i() * wo; degree]);

    // Cancel out gain change caused by inversion
    //     k_bs = k * real(prod(-z) / prod(-p));
    let k_bs = (zpk.k * prod( &negate(&zpk.z)) / prod( &negate(&zpk.p)) ).re;

    Ok( ZPKCoeffs {z: z_bs, p: p_bs, k: k_bs} )
}
// Lowpass
 fn lp2lp_zpk(zpk: &ZPKCoeffs, wo:f64) -> Result<ZPKCoeffs, Error>
 {
     let degree = relative_degree(zpk)?;

     // Scale all points radially from origin to shift cutoff frequency
     let z_lp:Vec<Complex<f64>> =
         zpk.z.iter()
             .map(|x| wo * x )
             .collect();
     let p_lp:Vec<Complex<f64>> =
         zpk.p.iter()
             .map(|x| wo * x )
             .collect();

     /* Each shifted pole decreases gain by wo, each shifted zero increases it.
        Cancel out the net change to keep overall gain the same. */
     let k_lp = zpk.k * wo.powf(degree as f64);


     Ok( ZPKCoeffs { z: z_lp, p: p_lp, k: k_lp} )
 }

//Highpass
fn lp2hp_zpk(zpk: &ZPKCoeffs, wo:f64) -> Result<ZPKCoeffs, Error>
{
    let degree = relative_degree(zpk)?;

    // Scale all points radially from origin to shift cutoff frequency
    let mut z_hp:Vec<Complex<f64>> =
        zpk.z.iter()
            .map(|x| wo / x )
            .collect();
    let p_hp:Vec<Complex<f64>> =
        zpk.p.iter()
            .map(|x| wo / x )
            .collect();

    // If lowpass had zeros at infinity, inverting moves them to origin.
    z_hp.append(&mut vec![Complex::<f64>::zero(); degree]);

    // Cancel out gain change caused by inversion
    let k_hp = zpk.k * (prod(&negate(&zpk.z)) / prod(&negate(&zpk.p)) ).re;

    Ok( ZPKCoeffs { z: z_hp, p: p_hp, k: k_hp} )
}

fn bilinear_zpk(zpk: &ZPKCoeffs, fs:f64) -> Result<ZPKCoeffs, Error>
{
    let one = Complex::<f64>::one();

    let degree = relative_degree(zpk)?;

    let z = zpk.z.clone();
    let p = zpk.p.clone();

    let fs2 = Complex::<f64>::new(2.0 * fs, 0.0);

    // Bilinear transform the poles and zeros
    
    // ndarray version: let mut z_z = (fs2 + &z) / (fs2 - &z);
    let mut z_z: Vec<Complex<f64>> = z.iter()
        .map(|z| (fs2 + z) / (fs2 - z))
        .collect();

    // ndarray version: let p_z = (fs2 + &p) / (fs2 - &p);
    let p_z:Vec<Complex<f64>> = p.iter()
        .map(|p| (fs2 + p) / (fs2 - p))
        .collect();

    // Any zeros that were at infinity get moved to the Nyquist frequency
    // ndarray version: z_z.append( Axis(0), (- Array::ones(degree)).view() )?;
    z_z.append(&mut vec![-one; degree]);

    // Compensate for gain change
    
    // ndarray version: let k_z = zpk.k * ( (fs2 - z).product() / (fs2 - p).product() ).re;
    let k_z = zpk.k * (
        z.iter()
            .map(|z| fs2 - z)
            .fold(one, |acc, x| acc * x) /
            p.iter()
                .map(|p| fs2 - p)
                .fold(one, |acc, x| acc * x))
        .re;

    Ok( ZPKCoeffs {z: z_z.to_vec(), p: p_z.to_vec(), k: k_z} )
}

///Return relative degree of transfer function from zeros and poles
fn relative_degree(zpk: &ZPKCoeffs) -> Result<usize, Error> {
 let degree = zpk.p.len() as i64 - zpk.z.len() as i64;
    if degree < 0 {
        return Err(InternalError("Improper transfer function. Must have at least as many poles as zeros.".to_string()));
    }

    Ok( degree as usize )
}

// Multiply all elements of the complex array and return the result.
// If array is empty return 1.
fn prod(list: &[Complex<f64>] ) -> Complex<f64> {

    let mut result = Complex::<f64>::one();

    if list.is_empty() {
        result
    } else {
        for c in list.iter() {
            result *= c;
        }

        result
    }
}

fn negate(list: &[Complex<f64>] ) -> Vec<Complex<f64>> {
    list.iter().map(|x| -x).collect()
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    use crate::vec_cplx;
    use crate::cplx;
    use num_traits::identities::Zero;
    use num_complex::Complex;
    use crate::filter_design::butterap;
    use crate::filter_design::FilterType::{BandStop, HighPass};
    use num_complex::ComplexFloat;
    use crate::util::Keys;

    #[test]
    fn test_butterap() {
        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![],
            p: vec_cplx![
            (-0.30901699437494745, 0.951_056_516_295_153_5),
            (-0.809_016_994_374_947_5, 0.587_785_252_292_473_1),
            ( -1.0, 0.0),
            (-0.809_016_994_374_947_5, -0.587_785_252_292_473_1),
            (-0.30901699437494745, -0.951_056_516_295_153_5)],
            k: 1.0,
        };

        let zpk_out = butterap(5).expect("butterap failed.");

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_butter__basic_digital_z_plane() {
        for N in 0..25 {
            let zpk = butter_internal(N, HighPass(0.01),  None).expect("Call to butter failed");

            assert_eq!(zpk.z.len(), N as usize);

            // All zeros exactly at DC
            for e in zpk.z.iter() {
                assert_eq!(*e, Complex::<f64>::one());
            }

            // No poles outside unit circle
            for e in zpk.p.iter() {
                assert!( e.abs() <= 1.0 );
            }
        }
    }

    #[test]
    fn test_butter_highpass_high_even_order() {
        let zpk_out = butter_internal(28, HighPass(0.43), None).expect("Call to butter failed");

        let zpk_expected = ZPKCoeffs {
            z: vec![Complex::<f64>::one(); 28],
            p: vec_cplx![
                (2.068257195514592e-01, 9.238294351481734e-01),
                (2.068257195514592e-01, -9.238294351481734e-01),
                (1.874933103892023e-01, 8.269455076775277e-01),
                (1.874933103892023e-01, -8.269455076775277e-01),
                (1.717435567330153e-01, 7.383078571194629e-01),
                (1.717435567330153e-01, -7.383078571194629e-01),
                (1.588266870755982e-01, 6.564623730651094e-01),
                (1.588266870755982e-01, -6.564623730651094e-01),
                (1.481881532502603e-01, 5.802343458081779e-01),
                (1.481881532502603e-01, -5.802343458081779e-01),
                (1.394122576319697e-01, 5.086609000582009e-01),
                (1.394122576319697e-01, -5.086609000582009e-01),
                (1.321840881809715e-01, 4.409411734716436e-01),
                (1.321840881809715e-01, -4.409411734716436e-01),
                (1.262633413354405e-01, 3.763990035551881e-01),
                (1.262633413354405e-01, -3.763990035551881e-01),
                (1.214660449478046e-01, 3.144545234797277e-01),
                (1.214660449478046e-01, -3.144545234797277e-01),
                (1.104868766650320e-01, 2.771505404367791e-02),
                (1.104868766650320e-01, -2.771505404367791e-02),
                (1.111768629525075e-01, 8.331369153155753e-02),
                (1.111768629525075e-01, -8.331369153155753e-02),
                (1.125740630842972e-01, 1.394219509611784e-01),
                (1.125740630842972e-01, -1.394219509611784e-01),
                (1.147138487992747e-01, 1.963932363793666e-01),
                (1.147138487992747e-01, -1.963932363793666e-01),
                (1.176516491045901e-01, 2.546021573417188e-01),
                (1.176516491045901e-01, - 2.546021573417188e-01)
            ],
            k: 1.446671081817286e-06,
        };

        zpk_out.assert_approx_equal_to_with_sort(&zpk_expected, 1E-7, &[Keys::Im], &[Keys::Im]);
    }

    #[test]
    fn test_butter_highpass_high_odd_order() {
        let zpk_out = butter_internal(27, HighPass(0.56), None).expect("Call to butter failed");

        let zpk_expected = ZPKCoeffs {
            z: vec![Complex::<f64>::one(); 27],
            p: vec_cplx![
                (-1.772572785680147e-01, 9.276431102995948e-01),
                (-1.772572785680147e-01, -9.276431102995948e-01),
                (-1.600766565322114e-01, 8.264026279893268e-01),
                (-1.600766565322114e-01, -8.264026279893268e-01),
                (-1.461948419016121e-01, 7.341841939120078e-01),
                (-1.461948419016121e-01, -7.341841939120078e-01),
                (-1.348975284762046e-01, 6.493235066053785e-01),
                (-1.348975284762046e-01, -6.493235066053785e-01),
                (-1.256628210712206e-01, 5.704921366889227e-01),
                (-1.256628210712206e-01, -5.704921366889227e-01),
                (-1.181038235962314e-01, 4.966120551231630e-01),
                (-1.181038235962314e-01, -4.966120551231630e-01),
                (-1.119304913239356e-01, 4.267938916403775e-01),
                (-1.119304913239356e-01, -4.267938916403775e-01),
                (-1.069237739782691e-01, 3.602914879527338e-01),
                (-1.069237739782691e-01, -3.602914879527338e-01),
                (-1.029178030691416e-01, 2.964677964142126e-01),
                (-1.029178030691416e-01, -2.964677964142126e-01),
                (-9.978747500816100e-02, 2.347687643085738e-01),
                (-9.978747500816100e-02, -2.347687643085738e-01),
                (-9.743974496324025e-02, 1.747028739092479e-01),
                (-9.743974496324025e-02, -1.747028739092479e-01),
                (-9.580754551625957e-02, 1.158246860771989e-01),
                (-9.580754551625957e-02, -1.158246860771989e-01),
                (-9.484562207782568e-02, 5.772118357151691e-02),
                (-9.484562207782568e-02, -5.772118357151691e-02),
                (-9.452783117928215e-02, 0.0)
            ],
            k: 9.585686688851069e-09,
        };

        zpk_out.assert_approx_equal_to_with_sort(&zpk_expected, 1E-8, &[Keys::Im], &[Keys::Im]);
    }

    #[test]
    fn test_butter_bandpass() {
        let zpk_out = butter_internal(8, BandPass(0.25, 0.33), None).expect("Call to butter failed");

        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![
                (1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 0.0),
                (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0), (-1.0, 0.0)
            ],
            p: vec_cplx![
                (4.979909925436156e-01, 8.367609424799387e-01),
                (4.979909925436156e-01, -8.367609424799387e-01),
                (4.913338722555539e-01, 7.866774509868817e-01),
                (4.913338722555539e-01, -7.866774509868817e-01),
                (5.035229361778706e-01, 7.401147376726750e-01),
                (5.035229361778706e-01, -7.401147376726750e-01),
                (5.307617160406101e-01, 7.029184459442954e-01),
                (5.307617160406101e-01, -7.029184459442954e-01),
                (5.680556159453138e-01, 6.788228792952775e-01),
                (5.680556159453138e-01, -6.788228792952775e-01),
                (6.100962560818854e-01, 6.693849403338664e-01),
                (6.100962560818854e-01, -6.693849403338664e-01),
                (6.904694312740631e-01, 6.930501690145245e-01),
                (6.904694312740631e-01, -6.930501690145245e-01),
                (6.521767004237027e-01, 6.744414640183752e-01),
                (6.521767004237027e-01, -6.744414640183752e-01)
            ],
            k: 3.398854055800844e-08,
        };

        zpk_out.assert_approx_equal_to_with_sort(&zpk_expected, 1E-13, &[Keys::Im], &[Keys::Im]);
    }

    #[test]
    fn test_butter_bandstop() {
        let zpk_out = butter_internal(7, BandStop(0.45, 0.56), None).expect("Call to butter failed");

        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![
                (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01),
              (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01),
              (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01),
              (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01),
              (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01),
              (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01),
              (-1.594474531383421e-02, 9.998728744679880e-01),
              (-1.594474531383421e-02, -9.998728744679880e-01)
            ],
            p: vec_cplx![
                (-1.766850742887729e-01, 9.466951258673900e-01),
                (-1.766850742887729e-01, -9.466951258673900e-01),
                (1.467897662432886e-01, 9.515917126462422e-01),
                (1.467897662432886e-01, -9.515917126462422e-01),
                (-1.370083529426906e-01, 8.880376681273993e-01),
                (-1.370083529426906e-01, -8.880376681273993e-01),
                (1.086774544701390e-01, 8.915240810704319e-01),
                (1.086774544701390e-01, -8.915240810704319e-01),
                (-7.982704457700891e-02, 8.506056315273435e-01),
                (-7.982704457700891e-02, -8.506056315273435e-01),
                (5.238812787110331e-02, 8.524011102699969e-01),
                (5.238812787110331e-02, -8.524011102699969e-01),
                (-1.357545000491310e-02, 8.382287744986582e-01),
                (-1.357545000491310e-02, -8.382287744986582e-01)
            ],
            k: 4.577122512960063e-01,
        };

        zpk_out.assert_approx_equal_to_with_sort(&zpk_expected, 1E-13, &[Keys::Im], &[Keys::Im]);
    }

    #[test]
    fn test_lp2bp_zpk_basic() {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(0.0, -2.0) , (0.0, 2.0)],
            p: vec_cplx![(-0.75, 0.0), (-0.5, -0.5), (-0.5, 0.5) ],
            k: 3.0
        };

        let a = cplx!( -225.0, -8.0).sqrt();
        let b = cplx!( -225.0, 8.0).sqrt();

        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![ (0.0,-25.0), (0.0, 25.0), (0.0, 9.0), (0.0, -9.0), (0.0,0.0) ],
            p: vec![
                Complex::<f64>::new(-3.0, -6.0 * f64::sqrt(6.0)),
                cplx!( -2.0, -2.0) + b,
                cplx!( -2.0, 2.0) + a,
                Complex::<f64>::new(-3.0, 6.0 * f64::sqrt(6.0)),
                cplx!( -2.0, -2.0) - b,
                cplx!( -2.0, 2.0) - a,
            ],
            k: 24.0,
        };

        let zpk_out = lp2bp_zpk(&zpk_input, 15.0, 8.0).expect("lp2bp_zpk() failed");

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_lp2bp_zpk2() {
        let zero = Complex::<f64>::zero();

        let zpk_expected = ZPKCoeffs {
            z: vec![zero, zero, zero, zero, zero, ],
            p: vec_cplx![
            (-0.0362762775801293, -0.1888301421919864),
            (-0.12367163989565325, -0.19209967216269005),
            (-0.28720279868217247, -0.24789295432066377),
            (-0.12367163989565325, 0.19209967216269005),
            (-0.0362762775801293, 0.1888301421919864),
            (-0.1412248136695468, 0.735122328561757),
            (-0.3410322500361953, 0.5297268111277003),
            (-0.28720279868217247, 0.24789295432066377),
            (-0.3410322500361953, -0.5297268111277003),
            (-0.1412248136695468, -0.735122328561757)
        ],
            k: 0.0625307037478668,
        };

        let zpk_input = ZPKCoeffs {
            z: vec_cplx![],
            p: vec_cplx![
                (-0.30901699437494745, 0.951_056_516_295_153_5),
                (-0.809_016_994_374_947_5, 0.587_785_252_292_473_1),
                ( -1.0, 0.0),
                (-0.809_016_994_374_947_5, -0.587_785_252_292_473_1),
                (-0.30901699437494745, -0.951_056_516_295_153_5)
            ],
            k: 1.0
        };

        let bw = 0.5744055973643449;
        let wo = 0.3793894626537474;

        let zpk_out = lp2bp_zpk(&zpk_input, wo, bw).expect("lp2bp_zpk() failed");

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_lp2bs_zpk_basic() {
        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(0.0, -2.0) , (0.0, 2.0)],
            p: vec_cplx![(-0.75, 0.0), (-0.5, -0.5), (-0.5, 0.5) ],
            k: 3.0,
        };

        let zpk_out = lp2bs_zpk(&zpk_input, 35.0, 12.0).expect("lp2bp_zpk() failed");

        let a = cplx!( 1234.0, 0.0).sqrt();
        let b = cplx!( 129.0, 0.0).sqrt();

        let i = Complex::<f64>::i();

        let zpk_expected = ZPKCoeffs {
            z: vec![
                3.0 * i + a * i,
                -3.0 * i - a * i,
                3.0 * i - a * i,
                -3.0 * i + a * i,
                35.0 * i,
                -35.0 * i,
            ],
            p: vec![
                -3.0 * i * b - 8.0,
                (-6.0 + 6.0 * i) + (-1225.0 - 72.0 * i).sqrt(),
                (-6.0 - 6.0 * i) + (-1225.0 + 72.0 * i).sqrt(),
                3.0 * i * b - 8.0,
                (-6.0 + 6.0 * i) - (-1225.0 - 72.0 * i).sqrt(),
                (-6.0 - 6.0 * i) - (-1225.0 + 72.0 * i).sqrt(),
            ],
            k: 32.0,
        };

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_lp2lp_zpk_basic1() {
        let i = Complex::<f64>::i();

        let zpk_input = ZPKCoeffs {
            z: vec_cplx![],
            p: vec![(-1.0+i)/f64::sqrt(2.0), (-1.0-i)/f64::sqrt(2.0) ],
            k: 1.0,
        };

        let zpk_out = lp2lp_zpk(&zpk_input, 5.0).expect("lp2bp_zpk() failed");

        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![],
            p: zpk_input.p.iter().map(|x| 5.0 * x).collect(),
            k: 25.0,
        };

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_lp2lp_zpk_pseudo_chebyshev() {
        let _i = Complex::<f64>::i();

        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(0.0, -2.0) , (0.0, 2.0)],
            p: vec_cplx![(-0.75, 0.0), (-0.5, -0.5), (-0.5, 0.5) ],
            k: 3.0,
        };

        let zpk_out = lp2lp_zpk(&zpk_input, 20.0).expect("lp2bp_zpk() failed");

        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![(0.0, -40.0), (0.0, 40.0)],
            p: vec_cplx![(-15.0, 0.0), (-10.0, -10.0), (-10.0, 10.0)],
            k: 60.0,
        };

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_lp2hp_zpk_basic1() {
        let i = Complex::<f64>::i();

        let zpk_input = ZPKCoeffs {
            z: vec_cplx![],
            p: vec![(-1.0+i)/f64::sqrt(2.0), (-1.0-i)/f64::sqrt(2.0) ],
            k: 1.0,
        };

        let zpk_out = lp2hp_zpk(&zpk_input, 5.0).expect("lp2bp_zpk() failed");

        let zpk_expected= ZPKCoeffs {
            z: vec_cplx![(0.0, 0.0), (0.0, 0.0)],
            p: zpk_input.p.iter().rev().map(|x| 5.0 * x).collect(),
            k: 1.0,
        };

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_lp2hp_zpk2() {
        let _i = Complex::<f64>::i();

        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(0.0, -2.0) , (0.0, 2.0)],
            p: vec_cplx![(-0.75, 0.0), (-0.5, -0.5), (-0.5, 0.5) ],
            k: 3.0,
        };

        let zpk_out = lp2hp_zpk(&zpk_input, 6.0).expect("lp2bp_zpk() failed");

        let zpk_expected = ZPKCoeffs {
            z: vec_cplx![(0.0, 3.0), (0.0, -3.0), (0.0, 0.0)],
            p: vec_cplx![(-8.0, 0.0), (-6.0, 6.0), (-6.0, -6.0)],
            k: 32.0,
        };

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }

    #[test]
    fn test_bilinear() {
        let _i = Complex::<f64>::i();

        let zpk_input = ZPKCoeffs {
            z: vec_cplx![(0.0, -2.0) , (0.0, 2.0)],
            p: vec_cplx![(-0.75, 0.0), (-0.5, -0.5), (-0.5, 0.5) ],
            k: 3.0,
        };

        let zpk_out = bilinear_zpk(&zpk_input, 10.0).expect("bilinear_zpk() failed");

        let zpk_expected = ZPKCoeffs {
            z: vec![
                cplx!( 20.0, -2.0) / cplx!( 20.0, 2.0),
                cplx!( 20.0, 2.0) / cplx!( 20.0, -2.0),
                cplx!( -1.0, 0.0),
            ],
            p: vec![
                cplx!(77.0/83.0, 0.0),
                cplx!(39.0/2.0, -1.0/2.0) / cplx!(41.0/2.0, 1.0/2.0 ),
                cplx!(39.0/2.0, 1.0/2.0) / cplx!(41.0/2.0, -1.0/2.0),
            ],
            k: 9696.0 / 69803.0,
        };

        zpk_out.assert_approx_equal_to(&zpk_expected, 1E-12);
    }
}