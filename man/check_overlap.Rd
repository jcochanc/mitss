\name{check_overlap}
\alias{check_overlap}
\title{Check overlap in propensity score distributions}
\usage{
check_overlap(e_hat, W, threshold = 0.3)
}
\arguments{
\item{e_hat}{Vector of estimated propensity scores}

\item{W}{Binary treatment indicator}

\item{threshold}{Maximum allowed Jensen-Shannon divergence (default: 0.3)}
}
\value{
List containing:
  \item{jsd}{Jensen-Shannon divergence}
  \item{adequate_overlap}{Boolean indicating if overlap is adequate}
  \item{e_min}{Minimum overlapping propensity score}
  \item{e_max}{Maximum overlapping propensity score}
  \item{trimmed_indices}{Indices of units in overlapping region}
}
\description{
Checks for adequate overlap between treatment and control groups
based on propensity score distributions
}
\examples{
n <- 1000
X <- matrix(rnorm(n), ncol = 1)
e <- plogis(X[,1])
W <- rbinom(n, 1, e)

overlap_check <- check_overlap(e, W)
print(overlap_check$jsd)
print(overlap_check$adequate_overlap)
}
