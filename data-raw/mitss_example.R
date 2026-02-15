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

usethis::use_data(mitss_example, overwrite = TRUE)