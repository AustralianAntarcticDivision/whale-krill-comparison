#' Compute Confidence Interval for a Log-Normal Mean
#'
#' @description
#' Given the mean and standard error (SE) of a variable assumed to follow a
#' log-normal distribution, this function returns the lower and upper confidence
#' interval bounds on the *original* (non-log) scale.
#'
#' @param mean Numeric. Mean of the variable on the original (non-log) scale.
#' @param se   Numeric. Standard error of the variable on the original scale.
#' @param conf.level Numeric. Confidence level (default = 0.95).
#'
#' @return A named data frame with lower and upper confidence limits.
#'
#' @details
#' The function internally converts mean and SE on the original scale to the
#' parameters of the corresponding log-normal distribution:
#' \deqn{
#'   \mu = \log\left(\frac{mean^2}{\sqrt{se^2 + mean^2}}\right),
#'   \sigma = \sqrt{\log\left(1 + \frac{se^2}{mean^2}\right)}
#' }
#' Confidence limits are calculated on the log scale and back-transformed
#' using the exponential function.
#'
#' @examples
#' lognormal_ci_orig(mean = 10, se = 2)               # 95% CI (default)
#' lognormal_ci_orig(mean = 10, se = 2, conf.level = 0.90)  # 90% CI
#'
lognormal_ci_orig <- function(mean, se, conf.level = 0.95) {
  # Check inputs
  #if (mean <= 0) stop("Mean must be positive for a log-normal distribution.")
  #if (se <= 0) stop("Standard error must be positive.")
  #if (conf.level <= 0 || conf.level >= 1)
  #  stop("Confidence level must be between 0 and 1.")
  
  # z-value for desired confidence interval
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  
  # Convert mean and se on original scale to parameters of lognormal
  sigma2 <- log(1 + (se^2 / mean^2))
  sigma <- sqrt(sigma2)
  mu <- log(mean) - 0.5 * sigma2
  
  # Compute lower and upper bounds on log scale
  lower_log <- mu - z * sigma
  upper_log <- mu + z * sigma
  
  # Transform back to original (non-log) scale
  lower <- exp(lower_log)
  upper <- exp(upper_log)
  
  # Return as named vector
  data.frame(lower = lower, upper = upper)
}
