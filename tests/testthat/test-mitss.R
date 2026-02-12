test_that("mitss works with simple data", {
  set.seed(123)
  n <- 200
  X <- rnorm(n)
  e <- plogis(0.5 * X)
  W <- rbinom(n, 1, e)
  Y <- 1.0 * W + 0.5 * X + rnorm(n)
  
  result <- mitss(Y, W, X, M = 5, seed = 456)
  
  expect_s3_class(result, "mitss")
  expect_true(!is.null(result$estimate))
  expect_true(!is.null(result$std_error))
  expect_true(result$ci_upper > result$ci_lower)
  expect_equal(length(result$imputed_datasets), 5)
})

test_that("mitss works with multiple covariates", {
  set.seed(789)
  n <- 300
  X <- matrix(rnorm(n * 3), ncol = 3)
  e <- plogis(0.3 * X[,1] - 0.2 * X[,2])
  W <- rbinom(n, 1, e)
  Y <- 0.5 * W + 0.3 * X[,1] + rnorm(n)
  
  result <- mitss(Y, W, X, M = 10, seed = 101)
  
  expect_s3_class(result, "mitss")
  expect_true(is.numeric(result$estimate))
  expect_true(result$std_error > 0)
})

test_that("mitss handles user-provided propensity scores", {
  set.seed(234)
  n <- 150
  X <- rnorm(n)
  e_true <- plogis(X)
  W <- rbinom(n, 1, e_true)
  Y <- W + X + rnorm(n)
  
  result <- mitss(Y, W, X, e_hat = e_true, M = 5, seed = 567)
  
  expect_equal(result$propensity_scores, e_true)
})

test_that("check_overlap detects poor overlap", {
  set.seed(345)
  n <- 200
  X <- rnorm(n)
  # Strong effect -> poor overlap
  e <- plogis(3 * X)
  W <- rbinom(n, 1, e)
  
  # Suppress expected warning about poor overlap
  expect_warning(
    overlap <- check_overlap(e, W, threshold = 0.3),
    "Limited overlap detected"
  )
  
  expect_true(is.numeric(overlap$jsd))
  expect_true(is.logical(overlap$adequate_overlap))
  expect_true(overlap$e_max > overlap$e_min)
})

test_that("trim_by_propensity removes correct units", {
  set.seed(456)
  n <- 200
  X <- rnorm(n)
  e <- plogis(2 * X)
  W <- rbinom(n, 1, e)
  Y <- W + X + rnorm(n)
  
  trimmed <- trim_by_propensity(Y, W, matrix(X, ncol=1), e, method = "minmax")
  
  expect_true(length(trimmed$Y) <= n)
  expect_equal(length(trimmed$Y), length(trimmed$W))
  expect_equal(nrow(trimmed$X), length(trimmed$Y))
  expect_true(trimmed$n_removed >= 0)
})

test_that("calculate_standardized_bias works", {
  set.seed(567)
  n <- 200
  X <- matrix(rnorm(n * 2), ncol = 2)
  colnames(X) <- c("X1", "X2")
  W <- rbinom(n, 1, 0.5)
  
  bias_df <- calculate_standardized_bias(X, W)
  
  expect_true(is.data.frame(bias_df))
  expect_equal(nrow(bias_df), 2)
  expect_true(all(c("variable", "std_bias") %in% names(bias_df)))
})

test_that("estimate_harmful_proportion works", {
  set.seed(678)
  n <- 200
  X <- rnorm(n)
  e <- plogis(0.5 * X)
  W <- rbinom(n, 1, e)
  Y <- 0.5 * W + 0.3 * X + rnorm(n)
  
  result <- mitss(Y, W, X, M = 5, seed = 789)
  harmful <- estimate_harmful_proportion(result)
  
  expect_true(is.numeric(harmful$prop_harmful))
  expect_true(harmful$prop_harmful >= 0 && harmful$prop_harmful <= 1)
  expect_true(harmful$ci_lower >= 0 && harmful$ci_upper <= 1)
  expect_true(harmful$ci_lower <= harmful$ci_upper)
})

test_that("print and summary methods work", {
  set.seed(890)
  n <- 150
  X <- rnorm(n)
  e <- plogis(0.5 * X)
  W <- rbinom(n, 1, e)
  Y <- W + X + rnorm(n)
  
  result <- mitss(Y, W, X, M = 5, seed = 901)
  
  expect_output(print(result), "MITSS")
  expect_output(summary(result), "Treatment Effect")
})

test_that("input validation works", {
  set.seed(111)
  n <- 100
  X <- rnorm(n)
  W <- rbinom(n, 1, 0.5)
  Y <- W + X + rnorm(n)
  
  # Mismatched lengths
  expect_error(mitss(Y, W[-1], X))
  
  # Non-binary treatment
  W_continuous <- runif(n)
  expect_error(mitss(Y, W_continuous, X))
  
  # Wrong X dimensions
  X_wrong <- matrix(rnorm((n-10) * 2), ncol = 2)
  expect_error(mitss(Y, W, X_wrong))
})
