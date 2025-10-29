#' Impute right-censored data via inverse-CDF from the conditional survival
#'
#' @param x Numeric vector of observed/censored values.
#' @param delta Numeric vector (same length as \code{x}): 1 = observed, 0 = right-censored.
#' @param theta Numeric parameter vector (see Details for each \code{dist}).
#' @param dist One of:
#' \itemize{
#'   \item \code{"normal"}: \code{theta = c(mean, variance)}
#'   \item \code{"poisson"}: \code{theta = lambda}
#'   \item \code{"weibull"}: \code{theta = c(rate, shape)}  (scale = 1/rate)
#'   \item \code{"gamma"}: \code{theta = c(shape, rate)}
#'   \item \code{"generalized gamma"}: \code{theta = c(alpha, beta, kappa)}
#'   \item \code{"lognormal"}: \code{theta = c(meanlog, varlog)}
#'   \item \code{"loglogistic"}: \code{theta = c(scale, shape)}
#' }
#' @return Numeric vector of imputations with observed entries unchanged.
#' @details
#' For each censored value \eqn{c}, draws \eqn{U \sim \mathrm{Unif}(F(c),1)} and returns
#' \eqn{F^{-1}(U)} of the chosen model (strict right-censoring, i.e., \eqn{T>c}).
#' @export
imputeRightCensor <- function(x, delta, theta, dist) {
  stopifnot(length(delta) == length(x))
  if (!all(delta %in% c(0, 1))) stop("delta must be 0/1.")
  if (!all(is.finite(x))) stop("x must be finite.")
  one_minus <- 1 - .Machine$double.eps

  t.imp <- x
  n <- length(x)

  if (dist == "normal") {
    mu <- theta[1]; sigma2 <- theta[2]
    if (!(sigma2 > 0)) stop("Normal requires variance > 0.")
    sd <- sqrt(sigma2)
    for (i in seq_len(n)) if (delta[i] == 0) {
      Fi <- pnorm(x[i], mean = mu, sd = sd)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      t.imp[i] <- qnorm(u, mean = mu, sd = sd)
    }

  } else if (dist == "poisson") {
    lambda <- theta[1]
    if (!(lambda > 0)) stop("Poisson requires lambda > 0.")
    for (i in seq_len(n)) if (delta[i] == 0) {
      Fi <- ppois(x[i], lambda = lambda)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      t.imp[i] <- qpois(u, lambda = lambda)
    }

  } else if (dist == "weibull") {
    rate <- theta[1]; shape <- theta[2]
    if (!(rate > 0 && shape > 0)) stop("Weibull requires rate > 0 and shape > 0.")
    scale <- 1 / rate
    for (i in seq_len(n)) if (delta[i] == 0) {
      Fi <- pweibull(x[i], shape = shape, scale = scale)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      t.imp[i] <- qweibull(u, shape = shape, scale = scale)
    }

  } else if (dist == "gamma") {
    shape <- theta[1]; rate <- theta[2]
    if (!(shape > 0 && rate > 0)) stop("Gamma requires shape > 0 and rate > 0.")
    for (i in seq_len(n)) if (delta[i] == 0) {
      Fi <- pgamma(x[i], shape = shape, rate = rate)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      t.imp[i] <- qgamma(u, shape = shape, rate = rate)
    }

  } else if (dist == "generalized gamma") {
    alpha <- theta[1]; beta <- theta[2]; kappa <- theta[3]
    if (!(alpha > 0 && beta > 0 && kappa > 0))
      stop("Generalized gamma requires alpha > 0, beta > 0, kappa > 0.")
    for (i in seq_len(n)) if (delta[i] == 0) {
      y0 <- (beta * x[i])^kappa
      Fi <- pgamma(y0, shape = alpha, rate = 1)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      # sample Y = Gamma^{-1}(u), then back-transform to T
      t.imp[i] <- .qquantile_ggamma(u, alpha = alpha, beta = beta, kappa = kappa)
    }

  } else if (dist == "lognormal") {
    meanlog <- theta[1]; varlog <- theta[2]
    if (!(varlog > 0)) stop("Lognormal requires varlog > 0.")
    sdlog <- sqrt(varlog)
    for (i in seq_len(n)) if (delta[i] == 0) {
      Fi <- plnorm(x[i], meanlog = meanlog, sdlog = sdlog)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      t.imp[i] <- qlnorm(u, meanlog = meanlog, sdlog = sdlog)
    }

  } else if (dist == "loglogistic") {
    scale <- theta[1]; shape <- theta[2]
    if (!(scale > 0 && shape > 0)) stop("Log-logistic requires scale > 0 and shape > 0.")
    for (i in seq_len(n)) if (delta[i] == 0) {
      # F(c) = (c/eta)^xi / (1 + (c/eta)^xi)
      a  <- (x[i] / scale)^shape
      Fi <- a / (1 + a)
      Fi <- min(Fi, one_minus)
      u  <- runif(1L, Fi, 1)
      t.imp[i] <- .qllogis(u, scale = scale, shape = shape)
    }

  } else {
    stop("Unsupported 'dist': ", dist,
         ". Supported: normal, poisson, weibull, gamma, generalized gamma, lognormal, loglogistic.")
  }

  t.imp
}
