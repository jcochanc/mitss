#!/usr/bin/env Rscript
# Demonstration of the mitss package

cat("MITSS Package Demonstration\n")
cat("============================\n\n")

# Source all R files
source("/mitss/R/mitss.R")
source("/mitss/R/methods.R")
source("/mitss/R/utilities.R")

# Example 1: Simple case with one covariate
cat("Example 1: Single Covariate\n")
cat("----------------------------\n")

set.seed(123)
n <- 500
X <- rnorm(n)

# Treatment assignment depends on X
e <- plogis(0.5 * X)
W <- rbinom(n, 1, e)

# Outcome: true treatment effect = 1.0
Y <- 1.0 * W + 0.5 * X + rnorm(n, sd = 0.5)

cat("True treatment effect: 1.0\n")
cat("Sample size:", n, "\n")
cat("Treatment proportion:", mean(W), "\n\n")

# Estimate treatment effect
result1 <- mitss(Y, W, matrix(X, ncol = 1), M = 10, seed = 456)

cat("MITSS Results:\n")
print(result1)
cat("\n")

# Example 2: Multiple covariates
cat("\nExample 2: Multiple Covariates\n")
cat("--------------------------------\n")

set.seed(789)
n <- 600
X <- matrix(rnorm(n * 3), ncol = 3)
colnames(X) <- paste0("X", 1:3)

# Treatment assignment
e <- plogis(0.3 * X[,1] - 0.4 * X[,2] + 0.2 * X[,3])
W <- rbinom(n, 1, e)

# Outcome: true treatment effect = 0.8
Y <- 2.0 + 0.8 * W + 0.3 * X[,1] + 0.2 * X[,2] - 0.1 * X[,3] + rnorm(n)

cat("True treatment effect: 0.8\n")
cat("Sample size:", n, "\n")
cat("Number of covariates:", ncol(X), "\n\n")

result2 <- mitss(Y, W, X, M = 15, seed = 101)

cat("MITSS Results:\n")
print(result2)
cat("\n")

# Check overlap
cat("\nChecking Propensity Score Overlap:\n")
cat("-----------------------------------\n")
overlap <- check_overlap(result2$propensity_scores, W)
cat("Jensen-Shannon Divergence:", round(overlap$jsd, 4), "\n")
cat("Adequate overlap (JSD < 0.3):", overlap$adequate_overlap, "\n")
cat("Propensity score range (treatment):", 
    round(range(result2$propensity_scores[W==1]), 3), "\n")
cat("Propensity score range (control):", 
    round(range(result2$propensity_scores[W==0]), 3), "\n\n")

# Covariate balance
cat("\nCovariate Balance Assessment:\n")
cat("------------------------------\n")
bias_df <- calculate_standardized_bias(X, W)
print(bias_df)
cat("\n")

# Non-linear estimand
cat("\nNon-linear Estimand: Proportion Harmed\n")
cat("---------------------------------------\n")
harmful <- estimate_harmful_proportion(result2)
cat("Proportion for whom treatment is not beneficial:\n")
cat("  Estimate:", round(harmful$prop_harmful, 3), "\n")
cat("  95% CI: [", round(harmful$ci_lower, 3), ",", 
    round(harmful$ci_upper, 3), "]\n\n")

# Comparison with naive estimator
naive_est <- mean(Y[W==1]) - mean(Y[W==0])
cat("\nComparison with Naive Estimator:\n")
cat("--------------------------------\n")
cat("Naive estimate (difference in means):", round(naive_est, 3), "\n")
cat("MITSS estimate:                      ", round(result2$estimate, 3), "\n")
cat("True effect:                         ", 0.8, "\n")
cat("MITSS bias:                          ", round(result2$estimate - 0.8, 3), "\n")
cat("Naive bias:                          ", round(naive_est - 0.8, 3), "\n\n")

cat("Demonstration complete!\n")
