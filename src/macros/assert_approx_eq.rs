/// Asserts that two expressions are equal to within a specified epsilon value.
///
/// An optional message can be specified with format arguments.
///
/// # Example
///
/// ```should_panic
/// use iir_filters::assert_approx_eq;
/// let temp = 65.0;
/// let expected_temp = 30.0;
/// let eps = 10.0;
///
/// assert_approx_eq!(temp, expected_temp, eps, "temp {} is incorrect! (Not within ϵ={} of {})", temp, eps, expected_temp);
/// ```
#[macro_export]
#[doc(hidden)]
macro_rules! assert_approx_eq {
    ($left:expr, $right:expr, $epsilon:expr $(,)?) => ({
        match (&$left, &$right, &$epsilon) {
            (left_val, right_val, epsilon) => {
                if ! ($crate::macros::util::eq_within_epsilon(*left_val, *right_val, *epsilon)) {
               // if ! ( (*left_val - *right_val).abs() <= *epsilon ) {
                    // The reborrows below are intentional. Without them, the stack slot for the
                    // borrow is initialized even before the values are compared, leading to a
                    // noticeable slow down.
                    panic!(r#"assertion failed: `(left == right)`
  left: `{:?}`,
 right: `{:?}`
   eps: `{:?}`"#, &*left_val, &*right_val, &*epsilon)
                }
            }
        }
    });
    ($left:expr, $right:expr, $epsilon:expr, $($arg:tt)+) => ({
        match (&($left), &($right), &($epsilon)) {
            (left_val, right_val, epsilon) => {
                //if !( (*left_val - *right_val).abs() <= *epsilon ) {
                 if ! ( $crate::macros::util::eq_within_epsilon(*left_val, *right_val, *epsilon)) {
                    // The reborrows below are intentional. Without them, the stack slot for the
                    // borrow is initialized even before the values are compared, leading to a
                    // noticeable slow down.
                    panic!(r#"assertion failed: `(left == right)`
  left: `{:?}`,
 right: `{:?}`
   eps: `{:?}`: {}"#, &*left_val, &*right_val, &*epsilon,
                           format_args!($($arg)+))
                }
            }
        }
    });
}

    #[cfg(test)]
    mod tests {
        #[test]
        fn compare_with_eps() {
            assert_approx_eq!(3f64, 4f64, 2f64);
        }

        #[test]
        #[should_panic]
        fn bad_compare_with_eps() {
            assert_approx_eq!(3f64, 4f64, 1e-3f64);
        }

        #[test]
        #[should_panic(expected = "assertion failed: `(left == right)`
  left: `3.0`,
 right: `4.0`
   eps: `0.001`: Temperature value incorrect.")]
        fn check_error_message() {
            assert_approx_eq!(3f64, 4f64, 1e-3f64, "Temperature value incorrect.");
        }

        // Make sure the value used for epsilon in the assert_eq
        // is the same as the value used in the error message.
        // PN: Not sure why this is relevant, but I updated this test to work with our implementation.
        #[test]
        #[should_panic(expected = "eps: `1.0")]
        fn should_evaluate_eps_only_once() {
            let mut count = 0_f64;

            // `count` will be 1.0 the first time the curly-braced block
            // is evaluated but 2.0 the second time.
            assert_approx_eq!(0_f64, 100_f64, {
            count += 1_f64;
            count
        });
    }
}
