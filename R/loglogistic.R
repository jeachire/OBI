#' Log-logistic density function
#'
#' Calculates the probability density function (PDF) for the log-logistic
#' distribution with scale parameter eta and shape parameter xi.
#'
#' @param t Vector of positive values where to evaluate the density.
#' @param eta Scale parameter (eta > 0).
#' @param xi Shape parameter (xi > 0).
#' @param log Logical; if TRUE, returns the logarithm of the density.
#'
#' @return Vector with density values (or log-density if log = TRUE).
#'
#' @export
#'
#' @examples
#' # Density at a point
#' dloglogistic(2, eta = 1.5, xi = 2)
dloglogistic <- function(t, eta, xi, log = FALSE) {
  # Parameter validation
  if (any(t <= 0)) stop("t must be greater than 0")
  if (eta <= 0) stop("eta must be greater than 0")
  if (xi <= 0) stop("xi must be greater than 0")

  # Density calculation
  term1 <- (xi/eta) * (t/eta)^(xi-1)
  term2 <- (1 + (t/eta)^xi)^2
  density <- term1 / term2

  # Apply log if requested
  if (log) {
    return(log(density))
  } else {
    return(density)
  }
}

#' Log-logistic cumulative distribution function
#'
#' Calculates the cumulative distribution function (CDF) for the log-logistic
#' distribution with scale parameter eta and shape parameter xi.
#'
#' @param t Vector of positive values where to evaluate the CDF.
#' @param eta Scale parameter (eta > 0).
#' @param xi Shape parameter (xi > 0).
#' @param lower.tail Logical; if TRUE (default), calculates P(T ≤ t),
#'                   if FALSE, calculates P(T > t) = 1 - F(t).
#' @param log.p Logical; if TRUE, returns the logarithm of the probability.
#'
#' @return Vector with cumulative distribution values
#'         (or log-probability if log.p = TRUE).
#'
#' @export
#'
#' @examples
#' # Cumulative probability at a point
#' ploglogistic(2, eta = 1.5, xi = 2)
ploglogistic <- function(t, eta, xi, lower.tail = TRUE, log.p = FALSE) {
  # Parameter validation
  if (any(t <= 0)) stop("t must be greater than 0")
  if (eta <= 0) stop("eta must be greater than 0")
  if (xi <= 0) stop("xi must be greater than 0")

  # CDF calculation
  cumulative <- 1 / (1 + (t/eta)^(-xi))

  # Adjust according to lower.tail
  if (!lower.tail) {
    cumulative <- 1 - cumulative
  }

  # Apply log if requested
  if (log.p) {
    return(log(cumulative))
  } else {
    return(cumulative)
  }
}
