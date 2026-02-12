\name{estimate_harmful_proportion}
\alias{estimate_harmful_proportion}
\title{Estimate non-linear estimands}
\usage{
estimate_harmful_proportion(mitss_result)
}
\arguments{
\item{mitss_result}{A mitss object from the mitss() function}
}
\value{
List containing:
  \item{prop_harmful}{Proportion for whom Y(1) <= Y(0)}
  \item{std_error}{Standard error of the proportion}
  \item{ci_lower}{Lower CI bound}
  \item{ci_upper}{Upper CI bound}
}
\description{
Estimates the proportion of units for whom treatment is harmful or not beneficial.
Uses the multiple imputation results from MITSS to properly quantify uncertainty.
}
\examples{
\dontrun{
# After running mitss()
result <- mitss(Y, W, X)
harmful <- estimate_harmful_proportion(result)
print(harmful$prop_harmful)
}
}
