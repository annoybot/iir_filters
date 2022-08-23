#[allow(clippy::expect_used)]
#[allow(clippy::panic)]
#[doc(hidden)]
pub(crate) mod test {
    use std::env;
    use std::str::FromStr;
    use crate::errors::Error;
    use crate::filter_design::ZPKCoeffs;
    use crate::sos::Sos;
    use crate::assert_approx_eq;
    use crate::util::{Keys, sort_lex};
    use crate::zpk2tf::BACoeffs;

    impl Sos {
        pub fn assert_approx_equal_to(&self, other: &Sos, eps: f64) {
            assert_eq!(self.num_sections(), other.num_sections(), "Number of sections are not equal.");

            for (i, (s1, s2)) in self.sections.iter().zip(&other.sections).enumerate() {
                assert_approx_eq!(s1.b0, s2.b0, eps, "b0 value at index {} not within ϵ", i);
                assert_approx_eq!(s1.b1, s2.b1, eps, "b1 value at index {} not within ϵ", i);
                assert_approx_eq!(s1.b2, s2.b2, eps, "b2 value at index {} not within ϵ", i);
                assert_approx_eq!(s1.a0, s2.a0, eps, "a0 value at index {} not within ϵ", i);
                assert_approx_eq!(s1.a1, s2.a1, eps, "a1 value at index {} not within ϵ", i);
                assert_approx_eq!(s1.a2, s2.a2, eps, "a2 value at index {} not within ϵ", i);
                
            }
        }
    }

    impl ZPKCoeffs {
        pub fn assert_approx_equal_to(&self, other: &ZPKCoeffs, eps: f64) {
            assert_eq!(self.z.len(), other.z.len(), "Number of zeroes are not equal.");
            assert_eq!(self.p.len(), other.p.len(), "Number of poles are not equal.");

            assert_approx_eq!(self.k, other.k, eps, "k value not within ϵ");

            for (i, (z1,z2)) in self.z.iter().zip(&other.z).enumerate() {
                assert_approx_eq!(*z1, *z2, eps, "z value at index {} not within ϵ", i);
            }

            for (i, (p1,p2)) in self.p.iter().zip(&other.p).enumerate() {
                assert_approx_eq!(*p1, *p2, eps, "p value at index {} not within ϵ", i);
            }
        }

        pub fn assert_approx_equal_to_with_sort(&self, other: &ZPKCoeffs, eps: f64, z_keys: &[Keys], p_keys: &[Keys]) {
            assert_eq!(self.z.len(), other.z.len(), "Number of zeroes are not equal.");
            assert_eq!(self.p.len(), other.p.len(), "Number of poles are not equal.");

            assert_approx_eq!(self.k, other.k, eps, "k value not within ϵ");

            let z1 = &mut self.z.clone();
            let z2 = &mut other.z.clone();

            sort_lex(z1, z_keys).expect("Sorting z1 failed.");
            sort_lex(z2, z_keys).expect("Sorting z2 failed.");

            for (i, (z1,z2)) in z1.iter().zip(z2).enumerate() {
                assert_approx_eq!(*z1, *z2, eps, "z value at index {} not within ϵ", i);
            }

            let p1 = &mut self.p.clone();
            let p2 = &mut other.p.clone();

            sort_lex(p1, p_keys).expect("Sorting p1 failed.");
            sort_lex(p2, p_keys).expect("Sorting p2 failed.");

            for (i, (p1,p2)) in p1.iter().zip(p2).enumerate() {
                assert_approx_eq!(*p1, *p2, eps, "p value at index {} not within ϵ", i);
            }
        }
    }

    impl BACoeffs {
        pub fn assert_approx_equal_to(&self, other: &BACoeffs, eps: f64) {
            assert_eq!(self.b.len(), other.b.len(), "Number of b's are not equal.");
            assert_eq!(self.a.len(), other.a.len(), "Number of a's are not equal.");

            for (i, (b1,b2)) in self.b.iter().zip(&other.b).enumerate() {
                assert_approx_eq!(*b1, *b2, eps, "b value at index {} not within ϵ", i);
            }

            for (i, (a1, a2)) in self.a.iter().zip(&other.a).enumerate() {
                assert_approx_eq!(*a1, *a2, eps, "a value at index {} not within ϵ", i);
            }
        }
    }

    // Routines related to reading test files.
    fn read_file<T: FromStr>(file_name: &str) -> Vec<Result<T, <T as FromStr>::Err>> {
        std::fs::read_to_string(file_name)
            .unwrap_or_else(|_| panic!("file {} not found!", file_name))
            .lines()
            .map(|x| x.parse())
            .collect()
    }

    #[derive(Debug, PartialEq)]
    struct DataPoint {
        t: f64,
        v: f64,
    }

    impl FromStr for DataPoint {
        type Err = Error;

        fn from_str(s: &str) -> Result<Self, Self::Err> {
            let x: Vec<&str> = s.split_whitespace().collect::<Vec<&str>>();
            let t = x[0].parse::<f64>()?;
            let v = x[1].parse::<f64>()?;

            Ok(DataPoint { t, v })
        }
    }
    //  Read one of our test files by name from the 'test_data' directory.
    //
    // Only the filename is required, this fn will figure out the full path.
    //Returns a tuple of two vecs, the first one contains time offsets, and the second contains the values.
    pub fn load_test_file(file_name: &str) -> Result<(Vec<f64>, Vec<f64>), Error> {
        let pathname = format!("{}/test_data/{}", env::var("CARGO_MANIFEST_DIR").expect("Env var CARGO_MANIFEST_DIR not found."), file_name );

        let parsed = read_file::<DataPoint>(&pathname);

        let mut t: Vec<f64> = vec![];
        let mut v: Vec<f64> = vec![];

        for e in parsed {
            let point = e?;

            t.push(point.t);
            v.push(point.v);
        }

        return Ok((t, v));
    }
}

// Test the code in the Readme when running `cargo test --doc`
#[cfg(doctest)]
mod test_readme {
    macro_rules! external_doc_test {
    ($x:expr) => {
        #[doc = $x]
        extern {}
    };
  }

    external_doc_test!(include_str!("../README.md"));
}