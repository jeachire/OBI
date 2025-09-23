#' @title Log-Likelihood Function Generator
#' @description Returns a function that computes the log-likelihood for a given
#' parametric lifetime distribution based on observed and censored data.
#'
#' @param dist A character string specifying the parametric distribution.
#' Supported options in the current implementation include:
#' \code{"normal"}, \code{"poisson"}, and \code{"weibull"}.
#'
#' @param x A numeric vector of data values (observed or censored).
#'
#' @param delta An optional numeric vector of the same length as \code{x},
#' indicating the censoring status of each observation:
#' \code{1} for observed values and \code{0} for right-censored values.
#' If \code{NULL}, all observations are assumed to be fully observed.
#'
#' @return A function of the parameter vector \code{theta}, which evaluates
#' the log-likelihood at the specified parameter values.
#'
#' @details
#' The generated log-likelihood function accounts for both observed and
#' right-censored data:
#' - For observed values (\code{delta = 1}), it uses the log-density.
#' - For censored values (\code{delta = 0}), it uses the log-survival function.
#'
#' Implemented distributions:
#' - \code{"normal"}: parameters = mean, sd
#' - \code{"poisson"}: parameter = lambda
#' - \code{"weibull"}: parameters = shape, scale
#'
#'
#' @examples
#' # Example with fully observed Normal data
#' loglik_fun <- getLogLikelihood("normal", x = rnorm(10, mean = 0, sd = 1))
#' loglik_fun(theta = c(0, 1))
#'
#' # Example with right-censored Weibull data
#' set.seed(123)
#' x <- rweibull(10, shape = 2, scale = 1)
#' delta <- rbinom(10, size = 1, prob = 0.7)  # 70% observed, 30% censored
#' loglik_fun <- getLogLikelihood("weibull", x = x, delta = delta)
#' loglik_fun(theta = c(2, 1))
#'
#' @export
getLogLikelihood <- function(dist, x, delta = NULL) {
  if (is.null(delta)) delta <- rep(1, length(x)) # all observed if not specified

  if (dist == "normal") {
    function(theta) {
      mu <- theta[1]; sigma <- theta[2]
      loglik_obs <- dnorm(x[delta == 1], mean = mu, sd = sigma, log = TRUE)
      loglik_cens <- pnorm(x[delta == 0], mean = mu, sd = sigma,
                           lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }
  } else if (dist == "poisson") {
    function(theta) {
      lambda <- theta[1]
      # Poisson is rarely censored, but support is included
      loglik_obs <- dpois(x[delta == 1], lambda, log = TRUE)
      loglik_cens <- ppois(x[delta == 0] - 1, lambda,
                           lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }
  } else if (dist == "weibull") {
    function(theta) {
      a <- theta[1]; sigma <- theta[2]
      loglik_obs <- dweibull(x[delta == 1], shape = a, scale = sigma, log = TRUE)
      loglik_cens <- pweibull(x[delta == 0], shape = a, scale = sigma,
                              lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }
  } else {
    stop("Unsupported distribution in getLogLikelihood().")
  }
}
