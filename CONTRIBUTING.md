# Thank you

If you are considering making a contribution to this project, thank you!

This document lays out some guidlines for contributing.  

The objective of this document isn't control, it's really to ensure you don't waste your time writing code which will not fit with the project's direction.

If you _do_ have a contribution you'd like to make that doesn't line up with these guidelines, please don't hesitate to open an issue so we can discuss it.

# Scope

This repo re-implements scipy's digital filters, so for the time being we will only accept contributions which derive froms scipy and which implement digital (not analog filters).

# Adding a new filter type

If you'd like to contribute a new filter type, these guidelines apply to you.

Suppose that you want to add Chebyshev type I digital filters. 


1. The first step is to find the scipy function which implements it, in this case: [cheby1()](https://github.com/scipy/scipy/blob/ee9985aafef8d3d99e19b2594767cf5fb8bfc3e6/scipy/signal/_filter_design.py#L3099) and implement it in [filter_design.rs](src/filter_design.rs).  Following the pattern for an existing algorithm, such as butterworth will give you an idea of how to proceed.

2. Almost more importantly than the code for desiging the filter, are the tests.  Locate any relevant [tests in the scypi source](https://github.com/scipy/scipy/blob/cbec0462607835cd38c2a03e80e7d23eba170aaa/scipy/signal/tests/test_filter_design.py#L2636) and implement all those which apply.  Since this crate doesn't implement the analog codepath, you can ignore those tests, but all other tests appying to digital filters for the new algorithm should be implemented.  Add a comment near your tests linking back to the python source.

3. Add a file based integration test.  These tests are defined in [filter.rs](src/filter.rs).  They are an end to end test which filters the file [test_data/bridge-strain-0-left-s-0-short.txt](test_data/bridge-strain-0-left-s-0-short.txt) and compares it the output the scipy routines.  This will involve writing a small python program to filter the same file and output the expected results in a file named after the unit test.

# Adding new dependencies

We would prefer not to add any additional dependencies to the crate.  We used to use `ndarray` which did make porting some scipy code easier, but we decided to follow an idiomatic Rust approach instead.

Again, if you have a reason to add a dependency which you think is persuasive, please open an issue so we can discuss it.