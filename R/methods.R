#' Print method for mitss objects
#'
#' @param x A mitss object
#' @param ... Additional arguments (not used)
#' @export
print.mitss <- function(x, ...) {
  cat("MITSS Causal Effect Estimate\n")
  cat("=============================\n\n")
  cat("Estimand:         ", x$estimand, "\n")
  cat("Point Estimate:   ", round(x$estimate, 4), "\n")
  cat("Standard Error:   ", round(x$std_error, 4), "\n")
  cat("95% CI:           [", round(x$ci_lower, 4), ", ", 
      round(x$ci_upper, 4), "]\n", sep = "")
  cat("\nNumber of Imputations: ", x$M, "\n")
  cat("Number of Subclasses:  ", x$n_subclasses, "\n")
  invisible(x)
}


#' Summary method for mitss objects
#'
#' @param object A mitss object
#' @param ... Additional arguments (not used)
#' @export
#' @method summary mitss
#' @importFrom stats median sd
summary.mitss <- function(object, ...) {
  cat("MITSS Causal Effect Estimation Summary\n")
  cat("=======================================\n\n")
  
  cat("Treatment Effect Estimate:\n")
  cat("  Point Estimate:  ", round(object$estimate, 4), "\n")
  cat("  Standard Error:  ", round(object$std_error, 4), "\n")
  cat("  95% CI:          [", round(object$ci_lower, 4), ", ", 
      round(object$ci_upper, 4), "]\n\n", sep = "")
  
  cat("Multiple Imputation Details:\n")
  cat("  Number of Imputations:  ", object$M, "\n")
  cat("  Range of Estimates:     [", round(min(object$individual_estimates), 4),
      ", ", round(max(object$individual_estimates), 4), "]\n", sep = "")
  cat("  SD of Estimates:        ", round(sd(object$individual_estimates), 4), "\n\n")
  
  cat("Propensity Score Summary:\n")
  cat("  Min:     ", round(min(object$propensity_scores), 4), "\n")
  cat("  Median:  ", round(median(object$propensity_scores), 4), "\n")
  cat("  Max:     ", round(max(object$propensity_scores), 4), "\n\n")
  
  cat("Subclass Information:\n")
  cat("  Number of Subclasses: ", object$n_subclasses, "\n")
  
  if (object$n_subclasses > 1) {
    subclass_table <- table(object$subclass_assignments)
    cat("  Units per Subclass:\n")
    for (i in 1:object$n_subclasses) {
      cat("    Subclass ", i, ": ", subclass_table[i], " units\n", sep = "")
    }
  }
  
  invisible(object)
}


#' Plot method for mitss objects
#'
#' @param x A mitss object
#' @param type Type of plot: "estimates" for distribution of MI estimates,
#'   "propensity" for propensity score distribution, or "balance" for covariate balance
#' @param ... Additional arguments passed to plot functions
#' @export
#' @importFrom graphics abline barplot hist legend rug
#' @importFrom stats density
plot.mitss <- function(x, type = c("estimates", "propensity", "balance"), ...) {
  type <- match.arg(type)
  
  if (type == "estimates") {
    # Histogram of individual imputation estimates
    hist(x$individual_estimates, 
         main = "Distribution of Treatment Effect Estimates Across Imputations",
         xlab = "Estimated Treatment Effect",
         col = "lightblue",
         border = "white",
         ...)
    abline(v = x$estimate, col = "red", lwd = 2, lty = 2)
    abline(v = x$ci_lower, col = "blue", lwd = 1, lty = 2)
    abline(v = x$ci_upper, col = "blue", lwd = 1, lty = 2)
    legend("topright", 
           legend = c("Point Estimate", "95% CI"),
           col = c("red", "blue"),
           lty = 2,
           lwd = c(2, 1))
    
  } else if (type == "propensity") {
    # Density plot of propensity scores
    plot(density(x$propensity_scores),
         main = "Distribution of Estimated Propensity Scores",
         xlab = "Propensity Score",
         ylab = "Density",
         ...)
    rug(x$propensity_scores)
    
  } else if (type == "balance") {
    # Simple balance plot showing subclass membership
    if (x$n_subclasses > 1) {
      barplot(table(x$subclass_assignments),
              main = "Distribution of Units Across Subclasses",
              xlab = "Subclass",
              ylab = "Number of Units",
              col = "lightgreen",
              ...)
    } else {
      message("No subclasses created - balance plot not available")
    }
  }
  
  invisible(x)
}


#' Extract coefficients from mitss object
#'
#' @param object A mitss object
#' @param ... Additional arguments (not used)
#' @export
#' @importFrom stats setNames
coef.mitss <- function(object, ...) {
  setNames(object$estimate, "Treatment Effect")
}


#' Extract confidence intervals from mitss object
#'
#' @param object A mitss object
#' @param parm Not used (for S3 consistency)
#' @param level Confidence level (default 0.95)
#' @param ... Additional arguments (not used)
#' @export
confint.mitss <- function(object, parm, level = 0.95, ...) {
  if (level != 0.95) {
    warning("Only 95% confidence intervals currently supported")
  }
  
  matrix(c(object$ci_lower, object$ci_upper),
         nrow = 1,
         dimnames = list("Treatment Effect", c("2.5 %", "97.5 %")))
}
