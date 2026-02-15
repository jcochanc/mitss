# mitss: Multiple Imputation using Two Subclassification Splines

<!-- badges: start -->
[![R-CMD-check](https://github.com/jcochanc/mitss/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jcochanc/mitss/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The `mitss` package implements the MITSS (Multiple Imputation using Two Subclassification Splines) method for estimating causal effects of binary treatments in observational studies. The method was developed by Gutman and Rubin (2015) and treats causal effect estimation as a missing data problem, using multiple imputation to properly quantify uncertainty.

## Key Features

- **Robust causal effect estimation** for binary treatments in unconfounded observational studies
- **Flexible response surface modeling** using splines with orthogonalized covariates
- **Proper uncertainty quantification** through multiple imputation
- **Diagnostic tools** for checking propensity score overlap and covariate balance
- **Non-linear estimands** such as the proportion of units harmed by treatment
- Works with both **scalar and multivariate covariates**

## Installation

You can install the package from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install mitss from GitHub
devtools::install_github("jcochanc/mitss", build_vignettes = TRUE)

# Or if you want to install without vignettes (faster)
devtools::install_github("jcochanc/mitss")
```

## Basic Usage

```r
library(mitss)

# Generate example data
set.seed(123)
n <- 500
X <- matrix(rnorm(n * 3), ncol = 3)
colnames(X) <- paste0("X", 1:3)

# Treatment assignment depends on covariates
e <- plogis(0.5 * X[,1] - 0.3 * X[,2])
W <- rbinom(n, 1, e)

# Outcome depends on treatment and covariates
Y <- 2 + 0.5 * W + 0.3 * X[,1] + 0.2 * X[,2] + rnorm(n)

# Estimate treatment effect
result <- mitss(Y, W, X, M = 25)

# View results
print(result)
summary(result)

# Visualize
plot(result, type = "estimates")
plot(result, type = "propensity")
```

## Checking Overlap

Before estimating treatment effects, it's important to check for adequate overlap in the propensity score distributions:

```r
# Estimate propensity scores (done automatically in mitss())
e_hat <- plogis(glm(W ~ X, family = binomial)$linear.predictors)

# Check overlap
overlap <- check_overlap(e_hat, W)
print(overlap$jsd)  # Jensen-Shannon divergence
print(overlap$adequate_overlap)  # TRUE if JSD < 0.3

# If overlap is poor, trim non-overlapping units
if (!overlap$adequate_overlap) {
  trimmed <- trim_by_propensity(Y, W, X, e_hat)
  result <- mitss(trimmed$Y, trimmed$W, trimmed$X, e_hat = trimmed$e_hat)
}
```

## Assessing Covariate Balance

```r
# Calculate standardized bias before and after subclassification
bias_before <- calculate_standardized_bias(X, W)
print(bias_before)

# After running mitss, check balance within subclasses
# (Can weight by propensity scores or use subclass indicators)
```

## Advanced Features

### Estimating Non-linear Estimands

```r
# Estimate proportion for whom treatment is harmful
harmful <- estimate_harmful_proportion(result)
cat("Proportion harmed:", harmful$prop_harmful, "\n")
cat("95% CI: [", harmful$ci_lower, ",", harmful$ci_upper, "]\n")
```

### Custom Propensity Score Models

```r
# Provide your own propensity score estimates
custom_e_hat <- # ... your propensity score estimates
result <- mitss(Y, W, X, e_hat = custom_e_hat)
```

### Adjusting Number of Subclasses

```r
# Use more subclasses (may improve fit but requires more data)
result <- mitss(Y, W, X, n_subclasses = 10)

# Or fewer subclasses (more robust with smaller samples)
result <- mitss(Y, W, X, n_subclasses = 4)
```

## Theoretical Background

The MITSS method views causal inference as a missing data problem. For each unit, we observe either Y(0) (potential outcome under control) or Y(1) (potential outcome under treatment), but not both. The method:

1. **Subclassifies** units based on estimated propensity scores to create balance
2. **Fits splines** across subclasses for flexible response surface modeling
3. **Orthogonalizes** covariates to reduce residual imbalance within subclasses
4. **Multiply imputes** missing potential outcomes
5. **Combines** results using Rubin's rules

This approach has been shown to have better coverage and efficiency compared to many standard methods like IPW, doubly robust estimation, and simple matching.

## Key Assumptions

- **Unconfoundedness**: Treatment assignment depends only on observed covariates
- **Overlap**: There is adequate overlap in the covariate distributions between treatment groups
- **No interference**: One unit's treatment doesn't affect another unit's outcome (SUTVA)

## Performance Considerations

- Works best when JSD < 0.3 (adequate overlap)
- Recommended: M = 25 imputations (higher for more precision)
- Minimum 3 units per treatment group per subclass
- Scales well to moderate-dimensional covariate spaces

## Citation

If you use this package, please cite:

Gutman, R. and Rubin, D.B. (2015). Estimation of causal effects of binary treatments in unconfounded studies. *Statistics in Medicine*, 34(26), 3381-3398. doi:10.1002/sim.6532

## References

- Gutman, R. and Rubin, D.B. (2015). Estimation of causal effects of binary treatments in unconfounded studies. *Statistics in Medicine*, 34(26), 3381-3398.
- Rubin, D.B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
- Imbens, G. and Rubin, D.B. (2015). *Causal Inference in Statistics, and in the Social and Biomedical Sciences*. Cambridge University Press.

## License

GPL (>= 3)

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.
