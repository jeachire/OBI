#' @title Impute Right-Censored Data
#' @description Imputes right-censored observations from a specified
#' parametric distribution, using random draws from the conditional
#' survival distribution.
#'
#' @param x A numeric vector of observed or censored survival times.
#'
#' @param delta A numeric vector of the same length as \code{x},
#' indicating the censoring status of each observation:
#' \code{1} for observed values and \code{0} for right-censored values.
#'
#' @param theta A numeric vector of distribution parameters.
#' The length and interpretation depend on \code{dist}:
#' - \code{"weibull"}: \code{c(shape, scale)}
#' - \code{"exponential"}: \code{c(rate)}
#'
#' @param dist A character string specifying the distribution.
#' Supported options are \code{"weibull"} and \code{"exponential"}.
#'
#' @return A numeric vector of the same length as \code{x}, where censored
#' observations (\code{delta = 0}) are replaced by random imputations
#' drawn from the conditional distribution, and observed values remain unchanged.
#'
#' @details
#' - For censored observations, the imputation is drawn from the
#' conditional distribution of the survival time given that it
#' exceeds the censoring time.
#' - For uncensored observations, the original values are preserved.
#'
#' The imputation is performed using the inversion method with random uniforms.
#'
#' @examples
#' set.seed(123)
#' # Example with Weibull distribution
#' x <- c(2, 5, 7, 3)
#' delta <- c(1, 0, 1, 0)  # two censored observations
#' theta <- c(shape = 2, scale = 3)
#' imputeRightCensor(x, delta, theta, dist = "weibull")
#'
#'
#' @export
imputeRightCensor <- function(x, delta, theta, dist = "weibull") {
  n <- length(x)
  t.imp <- x

  for (i in seq_len(n)) {
    if (delta[i] == 0) {
      if (dist == "weibull") {
        u <- runif(1, pweibull(x[i], shape = theta[1], scale = theta[2]), 1)
        t.imp[i] <- qweibull(u, shape = theta[1], scale = theta[2])
      } else if (dist == "exponential") {
        u <- runif(1, pexp(x[i], rate = theta[1]), 1)
        t.imp[i] <- qexp(u, rate = theta[1])
      } else {
        stop("Unsupported distribution in imputeRightCensor().")
      }
    }
  }
  return(t.imp)
}
