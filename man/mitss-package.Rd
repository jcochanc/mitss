\name{mitss-package}
\alias{mitss-package}
\docType{package}
\title{Multiple Imputation using Two Subclassification Splines for Causal Inference}
\description{
Implements the MITSS (Multiple Imputation using Two Subclassification 
Splines) method for estimating causal effects of binary treatments in 
unconfounded observational studies. Based on Gutman and Rubin (2015). 
The method treats causal effect estimation as a missing data problem 
and uses multiple imputation to estimate treatment effects with proper 
uncertainty quantification.
}
\details{
The main function is \code{\link{mitss}}, which estimates the average
treatment effect (ATE) using multiple imputation with spline-based
response surfaces.

Key features:
\itemize{
  \item Flexible response surface modeling with splines
  \item Automatic propensity score estimation
  \item Diagnostic tools for overlap and balance checking
  \item Support for non-linear estimands
  \item Proper uncertainty quantification via multiple imputation
}

Important functions:
\itemize{
  \item \code{\link{mitss}} - Main estimation function
  \item \code{\link{check_overlap}} - Check propensity score overlap
  \item \code{\link{trim_by_propensity}} - Remove non-overlapping units
  \item \code{\link{calculate_standardized_bias}} - Assess covariate balance
  \item \code{\link{estimate_harmful_proportion}} - Estimate proportion harmed
}
}
\author{
Roee Gutman, Donald B. Rubin (original method)

Maintainer: Your Name <your.email@example.com>
}
\references{
Gutman, R. and Rubin, D.B. (2015). Estimation of causal effects of binary 
treatments in unconfounded studies. Statistics in Medicine, 34(26), 3381-3398.
}
\keyword{package}
\seealso{
\code{\link{mitss}}
}
\examples{
# Basic usage
library(mitss)

# Generate example data
set.seed(123)
n <- 500
X <- matrix(rnorm(n * 3), ncol = 3)
colnames(X) <- paste0("X", 1:3)

e <- plogis(0.5 * X[,1] - 0.3 * X[,2])
W <- rbinom(n, 1, e)
Y <- 2 + 0.5 * W + 0.3 * X[,1] + 0.2 * X[,2] + rnorm(n)

# Estimate treatment effect
result <- mitss(Y, W, X, M = 10)
print(result)

# Check diagnostics
overlap <- check_overlap(result$propensity_scores, W)
balance <- calculate_standardized_bias(X, W)
}
