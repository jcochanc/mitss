\name{calculate_standardized_bias}
\alias{calculate_standardized_bias}
\title{Calculate standardized bias for covariates}
\usage{
calculate_standardized_bias(X, W, weights = NULL)
}
\arguments{
\item{X}{Covariate matrix}

\item{W}{Treatment indicator}

\item{weights}{Optional weights vector}
}
\value{
Data frame with standardized biases
}
\description{
Calculates standardized differences in covariate means between
treatment and control groups.
}
\examples{
set.seed(123)
n <- 500
X <- matrix(rnorm(n * 3), ncol = 3)
colnames(X) <- paste0("X", 1:3)
e <- plogis(0.5 * X[,1])
W <- rbinom(n, 1, e)

bias_df <- calculate_standardized_bias(X, W)
print(bias_df)
}
