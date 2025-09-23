#' @title Log-Prior Distribution Function Generator
#' @description Returns a function that computes the log-prior density
#' for the parameters of a specified distribution under the chosen prior.
#'
#' @param dist A character string specifying the parametric distribution.
#' Supported options in the current implementation include:
#' \code{"normal"}, \code{"poisson"}, and \code{"weibull"}.
#'
#' @param prior A character string or function specifying the prior distribution.
#' Built-in options include:
#' \code{"Jeffreys"} and \code{"uniform"} (case-sensitive).
#' Alternatively, the user may provide a custom function returning the
#' log-prior for a given parameter vector.
#'
#' @return A function of the parameter vector \code{theta},
#' which evaluates the log-prior density.
#'
#' @details
#' - \code{dist = "normal"} with \code{prior = "Jeffreys"}:
#'   Returns a constant prior (log-prior = 0), assuming known variance.
#' - \code{dist = "poisson"} with \code{prior = "Jeffreys"}:
#'   \deqn{\pi(\lambda) \propto \lambda^{-1/2}}
#'   so the log-prior is \eqn{-0.5 \log(\lambda)}.
#' - \code{dist = "weibull"} with \code{prior = "Jeffreys"}:
#'   \deqn{\pi(\alpha, \beta) \propto (\alpha \beta)^{-1}}
#'   so the log-prior is \eqn{-\log(\alpha) - \log(\beta)}.
#' - \code{prior = "uniform"}:
#'   Returns a flat prior, i.e., constant log-prior = 0.
#'
#' If a custom prior is passed as a function, the function is returned
#' directly without modification.
#'
#' @examples
#' # Jeffreys prior for Poisson distribution
#' logprior <- getPrior(dist = "poisson", prior = "Jeffreys")
#' logprior(theta = 3)
#'
#' # Jeffreys prior for Weibull distribution
#' logprior <- getPrior(dist = "weibull", prior = "Jeffreys")
#' logprior(theta = c(1.5, 2))
#'
#' # Uniform prior for Normal distribution
#' logprior <- getPrior(dist = "normal", prior = "uniform")
#' logprior(theta = c(0, 1))
#'
#' # Custom user-defined prior
#' myPrior <- function(theta) -sum(theta^2)  # log of a Gaussian-like prior
#' logprior <- getPrior(dist = "custom", prior = myPrior)
#' logprior(theta = c(1, 2))
#'
#' @export
getPrior <- function(dist, prior = "Jeffreys") {
  # If the user passes a custom prior function
  if (is.function(prior)) return(prior)

  # Built-in objective priors
  if (prior == "Jeffreys") {
    if (dist == "normal") return(function(theta) 0)
    if (dist == "poisson") return(function(theta) -0.5 * log(theta))
    if (dist == "weibull") return(function(theta) -log(theta[1]) - log(theta[2]))
  } else if (prior == "uniform") {
    return(function(theta) 0)
  } else {
    stop("Unknown prior: use 'Jeffreys', 'Reference', 'MDI', 'Uniform', or a custom function")
  }
}
