\name{trim_by_propensity}
\alias{trim_by_propensity}
\title{Trim units outside propensity score overlap region}
\usage{
trim_by_propensity(Y, W, X, e_hat, method = c("minmax", "crump"))
}
\arguments{
\item{Y}{Outcome vector}

\item{W}{Treatment vector}

\item{X}{Covariate matrix}

\item{e_hat}{Propensity scores}

\item{method}{Trimming method: "minmax" or "crump" (default: "minmax")}
}
\value{
List with trimmed data
}
\description{
Removes units that fall outside the propensity score overlap region
to improve the validity of causal effect estimates.
}
\examples{
set.seed(123)
n <- 500
X <- matrix(rnorm(n), ncol = 1)
e <- plogis(2 * X[,1])
W <- rbinom(n, 1, e)
Y <- W + X[,1] + rnorm(n)

trimmed <- trim_by_propensity(Y, W, X, e, method = "minmax")
cat("Trimmed", trimmed$n_removed, "units\\n")
}
