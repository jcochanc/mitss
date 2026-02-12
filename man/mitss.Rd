\name{mitss}
\alias{mitss}
\title{Multiple Imputation using Two Subclassification Splines (MITSS)}
\usage{
mitss(
  Y,
  W,
  X,
  e_hat = NULL,
  M = 25,
  n_subclasses = 6,
  min_per_subclass = 3,
  estimand = c("ATE", "finite_pop"),
  alpha = 0.05,
  seed = NULL
)
}
\arguments{
\item{Y}{Numeric vector of observed outcomes}

\item{W}{Binary vector of treatment assignments (0 = control, 1 = treatment)}

\item{X}{Matrix or data frame of covariates}

\item{e_hat}{Optional vector of estimated propensity scores. If NULL, will be estimated.}

\item{M}{Number of multiple imputations (default: 25)}

\item{n_subclasses}{Maximum number of subclasses (default: 6)}

\item{min_per_subclass}{Minimum units per treatment group in each subclass (default: 3)}

\item{estimand}{Type of estimand: "ATE" for average treatment effect, 
  "finite_pop" for finite population ATE (default: "ATE")}

\item{alpha}{Significance level for confidence intervals (default: 0.05)}

\item{seed}{Random seed for reproducibility (optional)}
}
\value{
A list containing:
  \item{estimate}{Point estimate of the treatment effect}
  \item{std_error}{Standard error of the estimate}
  \item{ci_lower}{Lower bound of confidence interval}
  \item{ci_upper}{Upper bound of confidence interval}
  \item{imputed_datasets}{List of M imputed datasets}
  \item{individual_estimates}{Vector of estimates from each imputation}
  \item{propensity_scores}{Estimated or provided propensity scores}
  \item{subclass_assignments}{Subclass assignments for each unit}
}
\description{
Estimates causal effects of binary treatments using multiple imputation
with spline-based response surface approximations.
}
\examples{
# Generate example data
set.seed(123)
n <- 500
X <- matrix(rnorm(n * 3), ncol = 3)
colnames(X) <- paste0("X", 1:3)

# Generate treatment assignment
e <- plogis(0.5 * X[,1] - 0.3 * X[,2])
W <- rbinom(n, 1, e)

# Generate outcomes
Y <- 2 + 0.5 * W + 0.3 * X[,1] + 0.2 * X[,2] + rnorm(n)

# Estimate treatment effect
result <- mitss(Y, W, X, M = 10)
print(result$estimate)
print(result$ci_lower)
print(result$ci_upper)
}
\references{
Gutman, R. and Rubin, D.B. (2015). Estimation of causal effects of binary 
treatments in unconfounded studies. Statistics in Medicine, 34(26), 3381-3398.
}
