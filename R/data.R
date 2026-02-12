#' Simulated observational study data
#'
#' A simulated dataset for demonstrating the MITSS method. This dataset
#' represents an observational study with a binary treatment, continuous
#' outcome, and multiple confounding covariates.
#'
#' @format A data frame with 800 observations and 6 variables:
#' \describe{
#'   \item{Y}{Continuous outcome variable}
#'   \item{W}{Binary treatment indicator (0 = control, 1 = treatment)}
#'   \item{X1}{Continuous covariate (standardized)}
#'   \item{X2}{Continuous covariate (standardized)}
#'   \item{X3}{Continuous covariate (standardized)}
#'   \item{propensity}{True propensity score (for simulation purposes)}
#' }
#'
#' @details
#' The data were generated using the following process:
#' \itemize{
#'   \item Covariates X1, X2, X3 are independently drawn from N(0,1)
#'   \item Propensity scores: e = logit^{-1}(0.3*X1 - 0.4*X2 + 0.2*X3)
#'   \item Treatment: W ~ Bernoulli(e)
#'   \item Outcome: Y = 2 + 0.8*W + 0.3*X1 + 0.2*X2 - 0.1*X3 + N(0,1)
#' }
#'
#' The true average treatment effect is 0.8.
#'
#' @examples
#' data(mitss_example)
#' head(mitss_example)
#'
#' # Estimate treatment effect
#' result <- mitss(
#'   Y = mitss_example$Y,
#'   W = mitss_example$W,
#'   X = as.matrix(mitss_example[, c("X1", "X2", "X3")]),
#'   M = 25
#' )
#' print(result)
#'
#' @source Simulated data based on Gutman and Rubin (2015)
"mitss_example"


# Generate the actual dataset
set.seed(14)
n <- 800
X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rnorm(n)

propensity <- plogis(0.3 * X1 - 0.4 * X2 + 0.2 * X3)
W <- rbinom(n, 1, propensity)

Y <- 2 + 0.8 * W + 0.3 * X1 + 0.2 * X2 - 0.1 * X3 + rnorm(n)

mitss_example <- data.frame(
  Y = Y,
  W = W,
  X1 = X1,
  X2 = X2,
  X3 = X3,
  propensity = propensity
)

# Save to data directory
# save(mitss_example, file = "data/mitss_example.rda", compress = "bzip2")
