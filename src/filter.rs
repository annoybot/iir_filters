//! Filter input using either the Direct Form 1 or Direct Form 2 Transposed representations.
//!
//! Adapted from crate [biquad](https://crates.io/crates/biquad).

use crate::sos::{Sos};

/// Internal states and coefficients of the Direct Form 1 form.
#[derive(Debug)]
pub struct DirectForm1 {
    y1: Vec<f64>,
    y2: Vec<f64>,
    x1: Vec<f64>,
    x2: Vec<f64>,
    sos: Sos,
}

/// Internal states and coefficients of the Direct Form 2 Transposed form.
#[derive(Debug)]
pub struct DirectForm2Transposed {
    s1: Vec<f64>,
    s2: Vec<f64>,
    sos: Sos,
}

/// The required functions of a filter implementation
pub trait Filter {
    /// Applies the filtering on the input.
    fn filter(&mut self, input: f64) -> f64;
}

impl DirectForm2Transposed {
    /// Creates a Direct Form 2 Transposed biquad from Sos filter coefficients
    pub fn new(sos: &Sos) -> Self {
        let num_sections = sos.sections.len();

        DirectForm2Transposed {
            s1: vec![0.0; num_sections],
            s2: vec![0.0; num_sections],
            sos: sos.clone(),
        }
    }
}

impl Filter for DirectForm2Transposed {
    /// Applies the filtering on the input, by cascading
    /// the second order sections.
    fn filter(&mut self, input: f64) -> f64 {
        let mut output = input;

        for (i, section) in self.sos.sections.iter().enumerate() {
            let prev_output = output;
            output = section.b0 * prev_output + self.s1[i] ;
            self.s1[i] =  section.b1 * prev_output - section.a1 * output + self.s2[i];
            self.s2[i] = section.b2 * prev_output - section.a2 * output;
        }

        output
    }
}

impl DirectForm1 {
    /// Creates a Direct Form 1 biquad from Sos filter coefficients
    pub fn new(sos: &Sos) -> Self {
        let num_sections = sos.sections.len();

        DirectForm1 {
            y1: vec![0.0;num_sections],
            y2: vec![0.0;num_sections],
            x1: vec![0.0;num_sections],
            x2: vec![0.0;num_sections],
            sos: sos.clone(),
        }
    }
}

impl Filter for DirectForm1 {
    /// Applies the filtering on the input, by cascading
    /// the second order sections.
    fn filter(&mut self, input: f64) -> f64 {
        let mut out = input;

        for (i, section) in self.sos.sections.iter().enumerate() {
            let prev_out = out;
            
            out = section.b0 * prev_out + section.b1 * self.x1[i] + section.b2 * self.x2[i]
                - section.a1 * self.y1[i]
                - section.a2 * self.y2[i];

            self.x2[i] = self.x1[i];
            self.x1[i] = prev_out;
            self.y2[i] = self.y1[i];
            self.y1[i] = out;
        }

        out
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use crate::filter_design::{butter, FilterType};
    use crate::sos::zpk2sos;
    use crate::filter::{DirectForm1, DirectForm2Transposed, Filter};
    use crate::test_util::test::load_test_file;
    use crate::assert_approx_eq;
    use stdext::function_name;

    enum FilterImpl {
        DF2T,
        DF1
    }

    #[rstest]
    #[case(FilterImpl::DF2T)]
    #[case(FilterImpl::DF1)]
    fn test_filter_butter_low_pass(#[case] filterImpl: FilterImpl) {
        let fn_name = function_name!();
        let results_filename = format!("{}_expected.txt", &fn_name[(fn_name.rfind(":").expect("Colon not found")+1)..]);

        // This test is based on the following on python example:
        //    sos = butter(5, [12], 'low', fs=81.0, output='sos')
        let order = 5;
        let cutoff = 12.0;
        let fs = 81.0;

        let zpk = butter(order, FilterType::LowPass(cutoff), fs).expect("butter failed.");
        let sos =zpk2sos(&zpk, None).expect("zpk2sos failed");

        let mut filter: Box<dyn Filter> = match filterImpl {
            FilterImpl::DF2T => { Box::new(DirectForm2Transposed::new(&sos)) }
            FilterImpl::DF1 => { Box::new(DirectForm1::new(&sos)) }
        };

        let (_, data)  = load_test_file("bridge-strain-0-left-s-0-short.txt").expect("Could not load file");
        let (_, expected_data) = load_test_file(&results_filename).expect("Could not load file");

        for (input, expected_output) in data.iter().zip(expected_data) {
            let output = filter.filter(*input);

            assert_approx_eq!(output, expected_output, 1E-12);
        }
    }

    #[rstest]
    #[case(FilterImpl::DF2T)]
    #[case(FilterImpl::DF1)]
    fn test_filter_butter_high_pass(#[case] filterImpl: FilterImpl) {
        let fn_name = function_name!();
        let results_filename = format!("{}_expected.txt", &fn_name[(fn_name.rfind(":").expect("Colon not found")+1)..]);

        // This test is based on the following on python example:
        //    sos = butter(5, [1], 'highpass', fs=81.0, output='sos')
        let order = 5;
        let cutoff = 1.0;
        let fs = 81.0;

        let zpk = butter(order, FilterType::HighPass(cutoff), fs).expect("butter failed.");
        let sos =zpk2sos(&zpk, None).expect("zpk2sos failed");

        let mut filter: Box<dyn Filter> = match filterImpl {
            FilterImpl::DF2T => { Box::new(DirectForm2Transposed::new(&sos)) }
            FilterImpl::DF1 => { Box::new(DirectForm1::new(&sos)) }
        };

        let (_, data)  = load_test_file("bridge-strain-0-left-s-0-short.txt").expect("Could not load file");
        let (_, expected_data) = load_test_file(&results_filename).expect("Could not load file");

        for (input, expected_output) in data.iter().zip(expected_data) {
            let output = filter.filter(*input);

            assert_approx_eq!(output, expected_output, 1E-12);
        }
    }

    #[rstest]
    #[case(FilterImpl::DF2T)]
    #[case(FilterImpl::DF1)]
    fn test_filter_butter_band_pass(#[case] filterImpl: FilterImpl) {
        let fn_name = function_name!();
        let results_filename = format!("{}_expected.txt", &fn_name[(fn_name.rfind(":").expect("Colon not found")+1)..]);

        // This test is based on the following on python example:
        //    sos = butter(5, [1, 10], 'bandpass', fs=81.0, output='sos')
        let order = 5;
        let cutoff_low = 1.0;
        let cutoff_hi= 10.0;
        let fs = 81.0;

        let zpk = butter(order,FilterType::BandPass(cutoff_low, cutoff_hi), fs).expect("butter failed.");
        let sos =zpk2sos(&zpk, None).expect("zpk2sos failed");

        let mut filter: Box<dyn Filter> = match filterImpl {
            FilterImpl::DF2T => { Box::new(DirectForm2Transposed::new(&sos)) }
            FilterImpl::DF1 => { Box::new(DirectForm1::new(&sos)) }
        };

        let (_, data)  = load_test_file("bridge-strain-0-left-s-0-short.txt").expect("Could not load file");
        let (_, expected_data) = load_test_file(&results_filename).expect("Could not load file");

        for (input, expected_output) in data.iter().zip(expected_data) {
            let output = filter.filter(*input);

            assert_approx_eq!(output, expected_output, 1E-12);
        }
    }

    #[rstest]
    #[case(FilterImpl::DF2T)]
    #[case(FilterImpl::DF1)]
    fn test_filter_butter_band_stop(#[case] filterImpl: FilterImpl) {
        let fn_name = function_name!();
        let results_filename = format!("{}_expected.txt", &fn_name[(fn_name.rfind(":").expect("Colon not found")+1)..]);

        // This test is based on the following on python example:
        //    sos = butter(5, [4, 5], 'bandstop', fs=81.0, output='sos')
        let order = 5;
        let cutoff_low = 4.0;
        let cutoff_hi= 5.0;
        let fs = 81.0;

        let zpk = butter(order,FilterType::BandStop(cutoff_low, cutoff_hi),fs).expect("butter failed.");
        let sos =zpk2sos(&zpk, None).expect("zpk2sos failed");

        let mut filter: Box<dyn Filter> = match filterImpl {
            FilterImpl::DF2T => { Box::new(DirectForm2Transposed::new(&sos)) }
            FilterImpl::DF1 => { Box::new(DirectForm1::new(&sos)) }
        };
        let (_, data)  = load_test_file("bridge-strain-0-left-s-0-short.txt").expect("Could not load file");
        let (_, expected_data) = load_test_file(&results_filename).expect("Could not load file");

        for (input, expected_output) in data.iter().zip(expected_data) {
            let output = filter.filter(*input);

            assert_approx_eq!(output, expected_output, 1E-12);
        }
    }
}

