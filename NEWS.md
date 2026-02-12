# mitss 0.1.0

## Initial Release

* Initial implementation of MITSS (Multiple Imputation using Two Subclassification Splines)
* Core function `mitss()` for estimating average treatment effects
* Support for both scalar and multivariate covariates
* Diagnostic functions:
  - `check_overlap()` for assessing propensity score overlap
  - `calculate_standardized_bias()` for covariate balance
  - `trim_by_propensity()` for removing non-overlapping units
* Non-linear estimand estimation with `estimate_harmful_proportion()`
* S3 methods: print, summary, plot, coef, confint
* Comprehensive documentation and vignettes
* Unit tests with testthat

## Reference

Based on:
Gutman, R. and Rubin, D.B. (2015). Estimation of causal effects of binary 
treatments in unconfounded studies. Statistics in Medicine, 34(26), 3381-3398.
