#' @title Log-Prior Distribution Function Generator
#' @description Returns a function that computes the log-prior density
#' for the parameters of a specified distribution under the chosen prior.
#'
#' @param dist Character string for the parametric distribution.
#' Supported options include: \code{"normal"}, \code{"poisson"}, \code{"gamma"},
#' \code{"generalized gamma"}, \code{"weibull"}, \code{"lognormal"}, \code{"loglogistic"}.
#' @param prior Character string for the prior: \code{"Jeffreys"}, \code{"reference"}, \code{"MDI"};
#' or a custom function returning a log-prior.
#'
#' @return A function of the parameter vector \code{theta} that evaluates the log-prior.
#'
#' @details
#' Parameterizations (your conventions):
#' \itemize{
#' \item \code{"normal"}: \code{theta = c(mean, variance)}
#' \item \code{"poisson"}: \code{theta = lambda}
#' \item \code{"weibull"}: \code{theta = c(rate, shape)}
#' \item \code{"gamma"}: \code{theta = c(shape, rate)}
#' \item \code{"generalized gamma"}: \code{theta = c(alpha, beta, kappa)}
#' \item \code{"lognormal"}: \code{theta = c(meanlog, varlog)}
#' \item \code{"loglogistic"}: \code{theta = c(scale, shape)}
#' }
#'
#' If a custom prior function is provided, it is returned as-is.
#' Some objective priors may lead to improper posteriors depending on the data.
#'
#' @examples
#' # Jeffreys prior for Weibull
#' logprior <- getPrior(dist = "weibull", prior = "Jeffreys")
#' logprior(theta = c(1.5, 2))
#'
#' # Custom prior
#' myPrior <- function(theta) -sum(theta^2)
#' logprior <- getPrior(dist = "weibull", prior = myPrior)
#' logprior(theta = c(1, 2))
#'
#' @export
getPrior <- function(dist, prior = "Jeffreys") {
  # Custom prior: return directly
  if (is.function(prior)) return(prior)

  # Normalize prior name (case-insensitive)
  prior_lc <- tolower(prior)
  valid_priors <- c("jeffreys", "reference", "mdi")
  if (!prior_lc %in% valid_priors) {
    stop("Unknown prior: use 'jeffreys', 'reference', 'mdi' or provide a custom function")
  }

  # Helper closures (for uniform domain checks)
  log_or_neg_inf <- function(cond, expr) if (cond) expr else -Inf

  if (prior_lc == "jeffreys") {

    if (dist == "normal") {
      return(function(theta) {
        sigma2 <- theta[2]
        log_or_neg_inf(sigma2 > 0, { -1.5*log(sigma2) })   # falta verificar pripiedad
      })
    }

    if (dist == "poisson") {
      return(function(theta) {
        lambda <- theta[1]
        log_or_neg_inf(lambda > 0, { -0.5 * log(lambda) })  # falta verificar propiedad
      })
    }

    if (dist == "gamma") {
      return(function(theta) {
        a <- theta[1]; b <- theta[2]
        log_or_neg_inf(a > 0 && b > 0, {
          0.5 * log(a * trigamma(a) - 1) - log(b)    # propia
        })
      })
    }

    if (dist == "weibull") {
      return(function(theta) {
        rate <- theta[1]; shape <- theta[2]
        log_or_neg_inf(rate > 0 && shape > 0, { -log(rate) - log(shape) })  # propia
      })
    }

    if (dist == "lognormal") {
      warning("Selected prior may yield an improper posterior for LogNormal (Jeffreys).", call. = FALSE)
      return(function(theta) {
        varlog <- theta[2]
        log_or_neg_inf(varlog > 0, { -log(varlog) })
      })
    }

    if (dist == "loglogistic") {
      return(function(theta) {
        eta <- theta[1]; xi <- theta[2]
        log_or_neg_inf(eta > 0 && xi > 0, { -log(eta) }) ## falta verificar propiedad
      })
    }

    if (dist == "generalized gamma") {
      warning("Selected prior may yield an improper posterior for GG distribution (Jeffreys).", call. = FALSE)
      return(function(theta) {
        a <- theta[1]; b <- theta[2]
        log_or_neg_inf(a > 0 && b > 0, {
          0.5 * log((a * trigamma(a))^2 - trigamma(a)- 1) - log(b)})
      })
    }
  }

  else if (prior_lc == "reference") {

    if (dist == "normal") {
      return(function(theta) {
        sigma2 <- theta[2]
        log_or_neg_inf(sigma2 > 0, { -log(sigma2)})  # falta verificar propiedad
      })
    }

    if (dist == "poisson") {
      return(function(theta) {
        lambda <- theta[1]
        log_or_neg_inf(lambda > 0, { -0.5 * log(lambda) })  # falta verificar propiedad
      })
    }

    if (dist == "gamma") {
      warning("Reference prior assuming beta is the parameter of interest")
      return(function(theta) {
        a <- theta[1]; b <- theta[2]
        log_or_neg_inf(a > 0 && b > 0, {
          0.5 * log(a * trigamma(a) - 1) - 0.5 * log(a) - log(b) # propia
          #  0.5 * log(trigamma(a)) - log(b)
        })
      })
    }

    if (dist == "weibull") {
      warning("Reference prior valid for both parameter orderings")
      return(function(theta) {
        rate <- theta[1]; shape <- theta[2]
        log_or_neg_inf(rate > 0 && shape > 0, { -log(rate) - log(shape) }) # proper
      })
    }

    if (dist == "lognormal") {
      return(function(theta) {
        varlog <- theta[2]
        log_or_neg_inf(varlog > 0, { -0.5 * log(varlog) })  # proper
      })
    }

    if (dist == "loglogistic") {
      return(function(theta) {
        scale <- theta[1]; shape <- theta[2]
        log_or_neg_inf(scale > 0 && shape > 0, { -log(scale) - log(shape) })  # falta verificar propiedad posterior
      })
    }

    if (dist == "generalized gamma") {
      warning("Reference prior derived for parameter ordering a > k > b; other orderings may lead to improper posteriors.", call. = FALSE)
      return(function(theta) {
        a <- theta[1]; b <- theta[2]; k <- theta[3]
        log_or_neg_inf(a > 0 && b > 0, {
          0.5 * log(a^2 +  (trigamma(a))^2 - trigamma(a)- 1) - log(b) - log(k) - 0.5 * log(a^2 * trigamma(a) + a - 1)})   # proper
      })
    }
  }

  else if (prior_lc == "mdi") {

    if (dist == "normal") {
      return(function(theta) {
        sigma2 <- theta[2]
        log_or_neg_inf(sigma2 > 0, { 0 })  # TODO: derive MDI normal
      })
    }

    if (dist == "poisson") {
      return(function(theta) {
        lambda <- theta[1]
        log_or_neg_inf(lambda > 0, { -0.5 * log(lambda) })  # TODO: derive MDI Poisson
      })
    }

    if (dist == "gamma") {
      return(function(theta) {
        a <- theta[1]; b <- theta[2]
        log_or_neg_inf(a > 0 && b > 0, {
          # TODO: verify MDI for gamma; current expression seems doubtful
          (a - 1) * digamma(a) / gamma(a) - a + log(b) - lgamma(a)
        })
      })
    }

    if (dist == "weibull") {
      warning("The posterior for Weibull with MDI prior may be improper")
      return(function(theta) {
        rate <- theta[1]; shape <- theta[2]
        log_or_neg_inf(rate > 0 && shape > 0, { log(rate) + log(shape) - rate/shape })  # TODO
      })
    }

    if (dist == "lognormal") {
      warning("The posterior for lognormal with MDI prior may be improper")
      return(function(theta) {
        varlog <- theta[2]
        log_or_neg_inf(varlog > 0, { -1.5 * log(varlog) })  # TODO: verify
      })
    }

    if (dist == "loglogistic") {
      warning("The posterior for loglogistic with MDI prior may be improper")
      return(function(theta) {
        scale <- theta[1]; shape <- theta[2]
        log_or_neg_inf(scale > 0 && shape > 0, { -log(scale) + log(shape) })  # TODO
      })
    }

    if (dist == "generalized gamma") {
      warning("MDI prior for generalized gamma not implemented; returning -Inf.")
      return(function(theta) -Inf)
    }
  }

  stop("Unsupported distribution: '", dist,
       "'. Supported: normal, poisson, gamma, generalized gamma, weibull, lognormal, loglogistic")
}
