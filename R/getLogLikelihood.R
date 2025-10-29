#' @title Log-Likelihood Function Generator
#' @description Returns a function that computes the log-likelihood for a given
#' parametric lifetime distribution based on observed and censored data.
#'
#' @param dist A character string specifying the parametric distribution.
#' Supported options in the current implementation include:
#' \code{"normal"}, \code{"poisson"}, \code{"gamma"},
#' \code{"weibull"}, \code{"generalized gamma"}, \code{"lognormal"},
#' and \code{"loglogistic"}.
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
#' \itemize{
#'   \item For observed values (\code{delta = 1}), it uses the log-density.
#'   \item For censored values (\code{delta = 0}), it uses the log-survival
#'   function, i.e., \eqn{\log S(x) = \log(1 - F(x))}.
#' }
#'
#' \strong{Censoring convention for Poisson:}
#' We use \emph{strict right-censoring}, so for censored counts the contribution is
#' \eqn{P(Y > x) = 1 - F(x)} (in R: \code{ppois(x, lambda, lower.tail = FALSE)}).
#'
#' \strong{Implemented distributions and parameter orders:}
#' \itemize{
#'   \item \code{"normal"}: parameters = \code{mean}, \code{variance}
#'   \item \code{"poisson"}: parameter = \code{lambda}
#'   \item \code{"weibull"}: parameters = \code{rate}, \code{shape}
#'   \item \code{"gamma"}: parameters = \code{shape}, \code{rate}
#'   \item \code{"generalized gamma"}: parameters = \code{alpha}, \code{beta}, \code{kappa}
#'   \item \code{"lognormal"}: parameters = \code{meanlog}, \code{variancelog}
#'   \item \code{"loglogistic"}: parameters = \code{scale}, \code{shape}
#' }
#'
#' \strong{Generalized gamma (Stacy) used here:}
#' With \eqn{\alpha>0}, \eqn{\beta>0}, \eqn{\kappa>0}, the density is
#' \deqn{
#' f(t \mid \alpha,\beta,\kappa) =
#' \frac{\kappa \beta^{\kappa \alpha}}{\Gamma(\alpha)}\,
#' t^{\kappa \alpha - 1}\exp\!\left[-(\beta t)^{\kappa}\right],\quad t>0.
#' }
#' The CDF is
#' \deqn{
#' F(t)=\frac{\gamma\!\left(\alpha,\ (\beta t)^{\kappa}\right)}{\Gamma(\alpha)}
#' \;=\; \mathrm{pgamma}\!\left((\beta t)^{\kappa};\ \text{shape}=\alpha,\ \text{rate}=1\right),
#' }
#' so we implement \code{dggamma()} and \code{pggamma()} accordingly and use
#' \code{pggamma(x, ..., lower.tail = FALSE, log.p = TRUE)} for censored terms.
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
#' # Example with generalized gamma (alpha, beta, kappa)
#' set.seed(1)
#' x <- rexp(20)  # solo para ejemplo; sustituye por tus datos
#' delta <- rbinom(20, 1, 0.8)
#' loglik_fun <- getLogLikelihood("generalized gamma", x = x, delta = delta)
#' # theta = c(alpha, beta, kappa)
#' loglik_fun(theta = c(2, 1, 1.5))
#'
#' @seealso dggamma, pggamma
#' @export

getLogLikelihood <- function(dist, x, delta = NULL) {
  if (is.null(delta)) delta <- rep(1L, length(x))
  stopifnot(is.numeric(x), is.numeric(delta), length(delta) == length(x))
  if (!all(delta %in% c(0,1))) stop("delta must be 0/1.")
  if (any(!is.finite(x))) stop("x must be finite.")

  if (dist == "normal") {
    function(theta) {
      mu <- theta[1]; sigma2 <- theta[2]
      if (!(sigma2 > 0)) return(-Inf)
      loglik_obs  <- dnorm(x[delta==1], mean = mu, sd = sqrt(sigma2), log = TRUE)
      loglik_cens <- pnorm(x[delta==0], mean = mu, sd = sqrt(sigma2),
                           lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else if (dist == "poisson") {
    function(theta) {
      lambda <- theta[1]
      if (!(lambda > 0)) return(-Inf)
      loglik_obs  <- dpois(x[delta == 1], lambda, log = TRUE)
      loglik_cens <- ppois(x[delta == 0], lambda, lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else if (dist == "weibull") {
    function(theta) {
      lambda <- theta[1]; kappa <- theta[2]
      if (!(lambda > 0 && kappa > 0)) return(-Inf)
      loglik_obs  <- dweibull(x[delta==1], scale = 1/lambda, shape = kappa, log = TRUE)
      loglik_cens <- pweibull(x[delta==0], scale = 1/lambda, shape = kappa,
                              lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else if (dist == "gamma") {
    function(theta) {
      alpha <- theta[1]; beta <- theta[2]
      if (!(alpha > 0 && beta > 0)) return(-Inf)
      loglik_obs  <- dgamma(x[delta==1], shape = alpha, rate = beta, log = TRUE)
      loglik_cens <- pgamma(x[delta==0], shape = alpha, rate = beta,
                            lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else if (dist == "generalized gamma") {
    function(theta) {
      alpha <- theta[1]; beta <- theta[2]; kappa <- theta[3]
      if (!(alpha > 0 && beta > 0 && kappa > 0)) return(-Inf)

      x_obs  <- x[delta == 1]
      x_cens <- x[delta == 0]

      loglik_obs  <- if (length(x_obs))  dggamma(x_obs,  alpha, beta, kappa, log = TRUE) else 0
      loglik_cens <- if (length(x_cens)) pggamma(x_cens, alpha, beta, kappa,
                                                 lower.tail = FALSE, log.p = TRUE) else 0
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else if (dist == "lognormal") {
    function(theta) {
      mu <- theta[1]; sigma2 <- theta[2]
      if (!(sigma2 > 0)) return(-Inf)
      loglik_obs  <- dlnorm(x[delta==1], meanlog = mu, sdlog = sqrt(sigma2), log = TRUE)
      loglik_cens <- plnorm(x[delta==0], meanlog = mu, sdlog = sqrt(sigma2),
                            lower.tail = FALSE, log.p = TRUE)
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else if (dist == "loglogistic") {
    function(theta) {
      eta <- theta[1]; xi <- theta[2]
      if (!(eta > 0 && xi > 0)) return(-Inf)
      x_obs  <- x[delta == 1]
      x_cens <- x[delta == 0]
      loglik_obs  <- if (length(x_obs))  dllogis_ll(x_obs,  scale = eta, shape = xi, log = TRUE) else 0
      loglik_cens <- if (length(x_cens)) pllogis_ll(x_cens, scale = eta, shape = xi,
                                                    lower.tail = FALSE, log.p = TRUE) else 0
      sum(loglik_obs) + sum(loglik_cens)
    }

  } else {
    stop("Unsupported distribution in getLogLikelihood().")
  }
}
