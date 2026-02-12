#' Check overlap in propensity score distributions
#'
#' Checks for adequate overlap between treatment and control groups
#' based on propensity score distributions
#'
#' @param e_hat Vector of estimated propensity scores
#' @param W Binary treatment indicator
#' @param threshold Maximum allowed Jensen-Shannon divergence (default: 0.3)
#'
#' @return List containing:
#'   \item{jsd}{Jensen-Shannon divergence}
#'   \item{adequate_overlap}{Boolean indicating if overlap is adequate}
#'   \item{e_min}{Minimum overlapping propensity score}
#'   \item{e_max}{Maximum overlapping propensity score}
#'   \item{trimmed_indices}{Indices of units in overlapping region}
#'
#' @export
#'
#' @examples
#' n <- 1000
#' X <- rnorm(n)
#' e <- plogis(X)
#' W <- rbinom(n, 1, e)
#' 
#' overlap_check <- check_overlap(e, W)
#' print(overlap_check$jsd)
#' print(overlap_check$adequate_overlap)
check_overlap <- function(e_hat, W, threshold = 0.3) {
  
  # Calculate Jensen-Shannon divergence
  jsd <- calculate_jsd(e_hat[W == 1], e_hat[W == 0])
  
  # Find overlapping region
  e_min <- max(min(e_hat[W == 0]), min(e_hat[W == 1]))
  e_max <- min(max(e_hat[W == 0]), max(e_hat[W == 1]))
  
  # Indices in overlapping region
  trimmed_indices <- which(e_hat >= e_min & e_hat <= e_max)
  
  adequate <- jsd < threshold
  
  if (!adequate) {
    warning(paste0("Limited overlap detected (JSD = ", round(jsd, 3), 
                   " > ", threshold, "). ",
                   "Consider trimming non-overlapping units."))
  }
  
  return(list(
    jsd = jsd,
    adequate_overlap = adequate,
    e_min = e_min,
    e_max = e_max,
    trimmed_indices = trimmed_indices,
    n_trimmed = length(e_hat) - length(trimmed_indices)
  ))
}


#' Calculate Jensen-Shannon divergence
#'
#' @param x1 Numeric vector (first distribution)
#' @param x2 Numeric vector (second distribution)
#' @param nbins Number of bins for histogram (default: 50)
#' @return Jensen-Shannon divergence
#' @keywords internal
calculate_jsd <- function(x1, x2, nbins = 50) {
  # Create common breaks
  all_x <- c(x1, x2)
  breaks <- seq(min(all_x), max(all_x), length.out = nbins + 1)
  
  # Get histograms
  h1 <- hist(x1, breaks = breaks, plot = FALSE)$counts
  h2 <- hist(x2, breaks = breaks, plot = FALSE)$counts
  
  # Convert to probabilities
  p <- h1 / sum(h1)
  q <- h2 / sum(h2)
  
  # Avoid log(0)
  p[p == 0] <- 1e-10
  q[q == 0] <- 1e-10
  
  # Calculate JSD
  lambda <- sum(h1) / (sum(h1) + sum(h2))
  m <- lambda * p + (1 - lambda) * q
  
  kl_pm <- sum(p * log(p / m))
  kl_qm <- sum(q * log(q / m))
  
  jsd <- lambda * kl_pm + (1 - lambda) * kl_qm
  
  return(jsd)
}


#' Trim units outside propensity score overlap region
#'
#' @param Y Outcome vector
#' @param W Treatment vector
#' @param X Covariate matrix
#' @param e_hat Propensity scores
#' @param method Trimming method: "minmax" or "crump" (default: "minmax")
#'
#' @return List with trimmed data
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' X <- rnorm(n)
#' e <- plogis(2 * X)
#' W <- rbinom(n, 1, e)
#' Y <- W + X + rnorm(n)
#' 
#' trimmed <- trim_by_propensity(Y, W, X, e, method = "minmax")
#' cat("Trimmed", trimmed$n_removed, "units\n")
trim_by_propensity <- function(Y, W, X, e_hat, method = c("minmax", "crump")) {
  method <- match.arg(method)
  
  if (method == "minmax") {
    # Keep units in [e_min, e_max]
    e_min <- max(min(e_hat[W == 0]), min(e_hat[W == 1]))
    e_max <- min(max(e_hat[W == 0]), max(e_hat[W == 1]))
    keep_idx <- which(e_hat >= e_min & e_hat <= e_max)
    
  } else if (method == "crump") {
    # Keep units with propensity scores in [0.1, 0.9]
    keep_idx <- which(e_hat >= 0.1 & e_hat <= 0.9)
  }
  
  n_removed <- length(Y) - length(keep_idx)
  
  if (n_removed > 0) {
    message(paste0("Trimmed ", n_removed, " units (", 
                   round(100 * n_removed / length(Y), 1), "%)"))
  }
  
  return(list(
    Y = Y[keep_idx],
    W = W[keep_idx],
    X = X[keep_idx, , drop = FALSE],
    e_hat = e_hat[keep_idx],
    keep_indices = keep_idx,
    n_removed = n_removed
  ))
}


#' Calculate standardized bias for covariates
#'
#' @param X Covariate matrix
#' @param W Treatment indicator
#' @param weights Optional weights vector
#'
#' @return Data frame with standardized biases
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 500
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' colnames(X) <- paste0("X", 1:3)
#' e <- plogis(0.5 * X[,1])
#' W <- rbinom(n, 1, e)
#' 
#' bias_df <- calculate_standardized_bias(X, W)
#' print(bias_df)
calculate_standardized_bias <- function(X, W, weights = NULL) {
  
  if (is.null(weights)) {
    weights <- rep(1, length(W))
  }
  
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  
  # Calculate weighted means and SDs
  mean_treat <- apply(X[W == 1, , drop = FALSE], 2, function(x) {
    weighted.mean(x, weights[W == 1])
  })
  
  mean_control <- apply(X[W == 0, , drop = FALSE], 2, function(x) {
    weighted.mean(x, weights[W == 0])
  })
  
  sd_treat <- apply(X[W == 1, , drop = FALSE], 2, function(x) {
    sqrt(weighted.mean((x - weighted.mean(x, weights[W == 1]))^2, 
                       weights[W == 1]))
  })
  
  sd_control <- apply(X[W == 0, , drop = FALSE], 2, function(x) {
    sqrt(weighted.mean((x - weighted.mean(x, weights[W == 0]))^2,
                       weights[W == 0]))
  })
  
  # Pooled SD
  sd_pooled <- sqrt((sd_treat^2 + sd_control^2) / 2)
  
  # Standardized bias
  std_bias <- (mean_treat - mean_control) / sd_pooled
  
  # Create output
  if (is.null(colnames(X))) {
    var_names <- paste0("X", 1:ncol(X))
  } else {
    var_names <- colnames(X)
  }
  
  result <- data.frame(
    variable = var_names,
    std_bias = std_bias,
    mean_treat = mean_treat,
    mean_control = mean_control,
    sd_pooled = sd_pooled
  )
  
  rownames(result) <- NULL
  
  return(result)
}


#' Estimate non-linear estimands
#'
#' Estimate the proportion of units for whom treatment is harmful or not beneficial
#'
#' @param mitss_result A mitss object from the mitss() function
#'
#' @return List containing:
#'   \item{prop_harmful}{Proportion for whom Y(1) <= Y(0)}
#'   \item{std_error}{Standard error of the proportion}
#'   \item{ci_lower}{Lower CI bound}
#'   \item{ci_upper}{Upper CI bound}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After running mitss()
#' result <- mitss(Y, W, X)
#' harmful <- estimate_harmful_proportion(result)
#' print(harmful$prop_harmful)
#' }
estimate_harmful_proportion <- function(mitss_result) {
  
  if (!inherits(mitss_result, "mitss")) {
    stop("Input must be a mitss object")
  }
  
  # Calculate proportion for each imputation
  props <- sapply(mitss_result$imputed_datasets, function(imp) {
    delta <- ifelse(imp$Y1_imputed - imp$Y0_imputed <= 0, 1, 0)
    mean(delta)
  })
  
  # Combine using MI rules
  M <- mitss_result$M
  prop_mean <- mean(props)
  
  # Within-imputation variance (binomial)
  W <- mean(props * (1 - props) / length(mitss_result$imputed_datasets[[1]]$Y0_imputed))
  
  # Between-imputation variance
  B <- var(props)
  
  # Total variance
  T_var <- W + (1 + 1/M) * B
  se <- sqrt(T_var)
  
  # CI
  df <- (M - 1) * (1 + W / ((1 + 1/M) * B))^2
  t_crit <- qt(0.975, df = df)
  
  return(list(
    prop_harmful = prop_mean,
    std_error = se,
    ci_lower = max(0, prop_mean - t_crit * se),
    ci_upper = min(1, prop_mean + t_crit * se),
    individual_props = props
  ))
}
