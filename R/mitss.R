#' Multiple Imputation using Two Subclassification Splines (MITSS)
#'
#' Estimates causal effects of binary treatments using multiple imputation
#' with spline-based response surface approximations.
#'
#' @param Y Numeric vector of observed outcomes
#' @param W Binary vector of treatment assignments (0 = control, 1 = treatment)
#' @param X Matrix or data frame of covariates
#' @param e_hat Optional vector of estimated propensity scores. If NULL, will be estimated.
#' @param M Number of multiple imputations (default: 25)
#' @param n_subclasses Maximum number of subclasses (default: 6)
#' @param min_per_subclass Minimum units per treatment group in each subclass (default: 3)
#' @param estimand Type of estimand: "ATE" for average treatment effect, 
#'   "finite_pop" for finite population ATE (default: "ATE")
#' @param alpha Significance level for confidence intervals (default: 0.05)
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A list containing:
#'   \item{estimate}{Point estimate of the treatment effect}
#'   \item{std_error}{Standard error of the estimate}
#'   \item{ci_lower}{Lower bound of confidence interval}
#'   \item{ci_upper}{Upper bound of confidence interval}
#'   \item{imputed_datasets}{List of M imputed datasets}
#'   \item{individual_estimates}{Vector of estimates from each imputation}
#'   \item{propensity_scores}{Estimated or provided propensity scores}
#'   \item{subclass_assignments}{Subclass assignments for each unit}
#'
#' @references
#' Gutman, R. and Rubin, D.B. (2015). Estimation of causal effects of binary 
#' treatments in unconfounded studies. Statistics in Medicine, 34(26), 3381-3398.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' colnames(X) <- paste0("X", 1:3)
#' 
#' # Generate treatment assignment
#' e <- plogis(0.5 * X[,1] - 0.3 * X[,2])
#' W <- rbinom(n, 1, e)
#' 
#' # Generate outcomes
#' Y <- 2 + 0.5 * W + 0.3 * X[,1] + 0.2 * X[,2] + rnorm(n)
#' 
#' # Estimate treatment effect
#' result <- mitss(Y, W, X, M = 10)
#' print(result$estimate)
#' print(result$ci_lower)
#' print(result$ci_upper)
#'
#' @export
mitss <- function(Y, W, X, e_hat = NULL, M = 25, n_subclasses = 6,
                  min_per_subclass = 3, estimand = c("ATE", "finite_pop"),
                  alpha = 0.05, seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Input validation
  estimand <- match.arg(estimand)
  n <- length(Y)
  
  if (length(W) != n) {
    stop("Y and W must have the same length")
  }
  
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  
  if (nrow(X) != n) {
    stop("X must have the same number of rows as length of Y")
  }
  
  if (!all(W %in% c(0, 1))) {
    stop("W must be binary (0 or 1)")
  }
  
  # Estimate propensity scores if not provided
  if (is.null(e_hat)) {
    e_hat <- estimate_propensity_score(W, X)
  }
  
  # Create subclasses
  subclass_info <- create_subclasses(e_hat, W, n_subclasses, min_per_subclass)
  subclasses <- subclass_info$subclasses
  n_subclasses_actual <- subclass_info$n_subclasses
  
  # Orthogonalize covariates
  X_ort <- orthogonalize_covariates(X, e_hat)
  
  # Prepare data
  data_list <- list(
    Y = Y,
    W = W,
    X = X,
    X_ort = X_ort,
    e_hat = e_hat,
    subclasses = subclasses
  )
  
  # Perform multiple imputation
  imputed_results <- replicate(M, {
    perform_single_imputation(data_list, n_subclasses_actual)
  }, simplify = FALSE)
  
  # Extract estimates and variances from each imputation
  gamma_m <- sapply(imputed_results, function(x) x$gamma)
  V_m <- sapply(imputed_results, function(x) x$V)
  
  # Combine results using Rubin's rules
  combined <- combine_mi_results(gamma_m, V_m, alpha)
  
  # Return results
  result <- list(
    estimate = combined$estimate,
    std_error = combined$std_error,
    ci_lower = combined$ci_lower,
    ci_upper = combined$ci_upper,
    imputed_datasets = imputed_results,
    individual_estimates = gamma_m,
    propensity_scores = e_hat,
    subclass_assignments = subclasses,
    estimand = estimand,
    M = M,
    n_subclasses = n_subclasses_actual
  )
  
  class(result) <- "mitss"
  return(result)
}


#' Estimate propensity scores
#'
#' @param W Binary treatment vector
#' @param X Covariate matrix
#' @return Vector of estimated propensity scores
#' @keywords internal
estimate_propensity_score <- function(W, X) {
  # Fit logistic regression
  X_df <- as.data.frame(X)
  formula_str <- paste("W ~", paste(names(X_df), collapse = " + "))
  
  # Add squared terms and interactions if X has multiple columns
  if (ncol(X) > 1) {
    # Add squared terms
    squared_terms <- paste0("I(", names(X_df), "^2)")
    formula_str <- paste(formula_str, "+", paste(squared_terms, collapse = " + "))
    
    # Add interactions for first few columns (to avoid over-parameterization)
    if (ncol(X) <= 5) {
      interaction_pairs <- combn(names(X_df), 2)
      interactions <- apply(interaction_pairs, 2, function(pair) {
        paste(pair[1], "*", pair[2])
      })
      formula_str <- paste(formula_str, "+", paste(interactions, collapse = " + "))
    }
  }
  
  formula_obj <- as.formula(formula_str)
  data_for_fit <- cbind(W = W, X_df)
  
  # Fit model with error handling
  fit <- tryCatch({
    glm(formula_obj, data = data_for_fit, family = binomial(link = "logit"))
  }, error = function(e) {
    # Fall back to simpler model if complex model fails
    simple_formula <- as.formula(paste("W ~", paste(names(X_df), collapse = " + ")))
    glm(simple_formula, data = data_for_fit, family = binomial(link = "logit"))
  })
  
  e_hat <- predict(fit, type = "response")
  
  # Bound away from 0 and 1
  e_hat <- pmax(pmin(e_hat, 0.99), 0.01)
  
  return(e_hat)
}


#' Create subclasses based on propensity scores
#'
#' @param e_hat Estimated propensity scores
#' @param W Treatment assignment
#' @param n_subclasses Maximum number of subclasses
#' @param min_per_subclass Minimum units per treatment group per subclass
#' @return List with subclass assignments and actual number of subclasses
#' @keywords internal
create_subclasses <- function(e_hat, W, n_subclasses, min_per_subclass) {
  n <- length(e_hat)
  
  # Start with quantile-based subclasses
  logit_e <- log(e_hat / (1 - e_hat))
  
  # Try to create n_subclasses, but may need fewer
  for (k in n_subclasses:2) {
    breaks <- quantile(logit_e, probs = seq(0, 1, length.out = k + 1))
    breaks[1] <- -Inf
    breaks[k + 1] <- Inf
    subclasses <- cut(logit_e, breaks = breaks, labels = FALSE)
    
    # Check if each subclass has minimum units per treatment group
    valid <- TRUE
    for (s in 1:k) {
      n_treat <- sum(subclasses == s & W == 1)
      n_control <- sum(subclasses == s & W == 0)
      if (n_treat < min_per_subclass | n_control < min_per_subclass) {
        valid <- FALSE
        break
      }
    }
    
    if (valid) {
      return(list(subclasses = subclasses, n_subclasses = k))
    }
  }
  
  # If we get here, even 2 subclasses doesn't work - return single subclass
  warning("Could not create subclasses with minimum units. Using single class.")
  return(list(subclasses = rep(1, n), n_subclasses = 1))
}


#' Orthogonalize covariates with respect to propensity score
#'
#' @param X Covariate matrix
#' @param e_hat Estimated propensity scores
#' @return Orthogonalized covariate matrix
#' @keywords internal
orthogonalize_covariates <- function(X, e_hat) {
  logit_e <- log(e_hat / (1 - e_hat))
  
  X_ort <- apply(X, 2, function(x) {
    residuals(lm(x ~ logit_e))
  })
  
  return(X_ort)
}


#' Perform single imputation
#'
#' @param data_list List containing Y, W, X, X_ort, e_hat, subclasses
#' @param n_subclasses Number of subclasses
#' @return List with imputed values and estimate
#' @keywords internal
perform_single_imputation <- function(data_list, n_subclasses) {
  Y <- data_list$Y
  W <- data_list$W
  X <- data_list$X
  X_ort <- data_list$X_ort
  e_hat <- data_list$e_hat
  subclasses <- data_list$subclasses
  
  n <- length(Y)
  logit_e <- log(e_hat / (1 - e_hat))
  
  # Fit response surfaces for each treatment group
  Y0_imputed <- rep(NA, n)
  Y1_imputed <- rep(NA, n)
  
  # For treatment group (W=1)
  idx_1 <- which(W == 1)
  if (length(idx_1) > 0) {
    fit_1 <- fit_spline_model(Y[idx_1], logit_e[idx_1], X_ort[idx_1, , drop = FALSE],
                              subclasses[idx_1], n_subclasses)
    
    # Observed values
    Y1_imputed[idx_1] <- Y[idx_1]
    
    # Impute missing Y(1) for control group
    idx_0 <- which(W == 0)
    if (length(idx_0) > 0) {
      Y1_imputed[idx_0] <- predict_from_spline(fit_1, logit_e[idx_0], 
                                               X_ort[idx_0, , drop = FALSE],
                                               subclasses[idx_0])
    }
  }
  
  # For control group (W=0)
  idx_0 <- which(W == 0)
  if (length(idx_0) > 0) {
    fit_0 <- fit_spline_model(Y[idx_0], logit_e[idx_0], X_ort[idx_0, , drop = FALSE],
                              subclasses[idx_0], n_subclasses)
    
    # Observed values
    Y0_imputed[idx_0] <- Y[idx_0]
    
    # Impute missing Y(0) for treatment group
    idx_1 <- which(W == 1)
    if (length(idx_1) > 0) {
      Y0_imputed[idx_1] <- predict_from_spline(fit_0, logit_e[idx_1],
                                               X_ort[idx_1, , drop = FALSE],
                                               subclasses[idx_1])
    }
  }
  
  # Calculate treatment effect
  gamma <- mean(Y1_imputed - Y0_imputed)
  
  # For super-population, V is estimated from imputation variance
  # For finite population, V = 0
  V <- 0
  
  return(list(
    Y0_imputed = Y0_imputed,
    Y1_imputed = Y1_imputed,
    gamma = gamma,
    V = V
  ))
}


#' Fit spline model with orthogonalized covariates
#'
#' @param Y Outcome values
#' @param logit_e Logit of propensity scores
#' @param X_ort Orthogonalized covariates
#' @param subclasses Subclass assignments
#' @param n_subclasses Number of subclasses
#' @return Fitted model object
#' @keywords internal
fit_spline_model <- function(Y, logit_e, X_ort, subclasses, n_subclasses) {
  # Create knots at subclass boundaries
  if (n_subclasses > 1) {
    knots <- numeric(n_subclasses - 1)
    for (k in 1:(n_subclasses - 1)) {
      # Knot at boundary between subclass k and k+1
      vals_k <- logit_e[subclasses == k]
      vals_kp1 <- logit_e[subclasses == (k + 1)]
      if (length(vals_k) > 0 && length(vals_kp1) > 0) {
        knots[k] <- mean(c(max(vals_k), min(vals_kp1)))
      } else {
        knots[k] <- quantile(logit_e, probs = k / n_subclasses)
      }
    }
  } else {
    knots <- NULL
  }
  
  # Create spline basis
  if (!is.null(knots) && length(knots) > 0) {
    spline_basis <- splines::bs(logit_e, knots = knots, degree = 3)
  } else {
    spline_basis <- splines::bs(logit_e, df = 4, degree = 3)
  }
  
  # Prepare design matrix
  if (ncol(X_ort) > 0) {
    design_matrix <- cbind(spline_basis, X_ort)
  } else {
    design_matrix <- spline_basis
  }
  
  # Fit linear model
  fit <- lm(Y ~ design_matrix)
  
  # Store information for prediction
  fit$knots <- knots
  fit$logit_e_range <- range(logit_e)
  fit$X_ort_included <- ncol(X_ort) > 0
  
  return(fit)
}


#' Predict from fitted spline model
#'
#' @param fit Fitted model object
#' @param logit_e_new Logit propensity scores for prediction
#' @param X_ort_new Orthogonalized covariates for prediction
#' @param subclasses_new Subclass assignments for prediction
#' @return Predicted values with random draws from posterior
#' @keywords internal
predict_from_spline <- function(fit, logit_e_new, X_ort_new, subclasses_new) {
  # Create spline basis for new data
  if (!is.null(fit$knots) && length(fit$knots) > 0) {
    spline_basis_new <- splines::bs(logit_e_new, knots = fit$knots, degree = 3)
  } else {
    # Use same df as in fitting
    spline_basis_new <- splines::bs(logit_e_new, df = 4, degree = 3)
  }
  
  # Prepare design matrix
  if (fit$X_ort_included && ncol(X_ort_new) > 0) {
    design_matrix_new <- cbind(spline_basis_new, X_ort_new)
  } else {
    design_matrix_new <- spline_basis_new
  }
  
  # Get point predictions
  newdata <- data.frame(design_matrix = I(design_matrix_new))
  point_pred <- predict(fit, newdata = newdata)
  
  # Draw from posterior predictive distribution
  sigma_hat <- summary(fit)$sigma
  n_new <- length(logit_e_new)
  
  # Add random variation
  predictions <- point_pred + rnorm(n_new, mean = 0, sd = sigma_hat)
  
  return(predictions)
}


#' Combine multiple imputation results using Rubin's rules
#'
#' @param gamma_m Vector of estimates from each imputation
#' @param V_m Vector of variances from each imputation
#' @param alpha Significance level
#' @return List with combined estimate, standard error, and CI
#' @keywords internal
combine_mi_results <- function(gamma_m, V_m, alpha) {
  M <- length(gamma_m)
  
  # Point estimate: average across imputations
  gamma_bar <- mean(gamma_m)
  
  # Within-imputation variance
  W <- mean(V_m)
  
  # Between-imputation variance
  B <- var(gamma_m)
  
  # Total variance
  T_var <- W + (1 + 1/M) * B
  
  # Standard error
  se <- sqrt(T_var)
  
  # Degrees of freedom (Barnard-Rubin adjustment)
  df <- (M - 1) * (1 + W / ((1 + 1/M) * B))^2
  
  # Confidence interval
  t_crit <- qt(1 - alpha/2, df = df)
  ci_lower <- gamma_bar - t_crit * se
  ci_upper <- gamma_bar + t_crit * se
  
  return(list(
    estimate = gamma_bar,
    std_error = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    df = df,
    within_var = W,
    between_var = B
  ))
}
