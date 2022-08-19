//! A Rust re-implementation of some of scipy's digital filters.
//!
//! # Example
//!
//! ```rust
//! use iir_filters::filter_design::FilterType;
//! use iir_filters::filter_design::butter;
//! use iir_filters::sos::zpk2sos;
//! use iir_filters::filter::DirectForm2Transposed;
//! use iir_filters::filter::Filter;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let order = 5;
//!     let cutoff_low = 1.0;
//!     let cutoff_hi= 10.0;
//!     let fs = 81.0;
//!
//!     let zpk = butter(order, FilterType::BandPass(cutoff_low, cutoff_hi),fs)?;
//!     let sos = zpk2sos(&zpk, None)?;
//!
//!     let mut dft2 = DirectForm2Transposed::new(&sos);
//!
//!     let data:Vec<f64>  = vec![1.0, 2.0, 3.0];
//!     let mut output:Vec<f64> = vec![];
//!
//!     for x in data.iter() {
//!         output.push( dft2.filter(*x) );
//!     }
//!
//!     return Ok( () );
//! }
//! ```
//! # Notes
//!
//! See: [scipy.signal: butter()](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html)
//!
//! ⚠️ For now it only implements Butterworth filters, because that's all I'm interested in.
//!

#![ cfg_attr( nightly, feature(doc_auto_cfg) ) ]

#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(clippy::needless_return)]
#![deny(unsafe_code)]
#![deny(unused_must_use)]
#![deny(clippy::panic)]
#![deny(clippy::expect_used)]
#![deny(clippy::unwrap_used)]
#![warn(missing_docs)]

pub mod filter_design;
pub mod errors;

#[doc(hidden)]
pub mod macros;

mod util;
pub mod sos;
mod cplxreal;
mod zpk2tf;
pub mod filter;
mod test_util;


