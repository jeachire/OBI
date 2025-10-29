#' @title Log-logistic density (scale–shape parameterization)
#' @name dllogis_ll
#' @description Computes the density of the log-logistic distribution for
#' \eqn{x>0} under the parameterization with \code{scale = \eqn{\eta} > 0}
#' and \code{shape = \eqn{\xi} > 0}.
#'
#' @details
#' For \eqn{x>0}, the density is
#' \deqn{
#' f(x;\eta,\xi) \;=\; \frac{\xi}{\eta}\,
#' \left(\frac{x}{\eta}\right)^{\xi-1}\,
#' \left[1+\left(\frac{x}{\eta}\right)^{\xi}\right]^{-2}.
#' }
#' For \eqn{x \le 0}, \code{dllogis_ll} returns 0 (or \code{-Inf} if \code{log=TRUE}).
#'
#' @param x Numeric vector of quantiles (\eqn{x}).
#' @param scale Positive scale parameter \eqn{\eta > 0}.
#' @param shape Positive shape parameter \eqn{\xi > 0}.
#' @param log Logical; if \code{TRUE}, return the log-density.
#'
#' @return A numeric vector of the same length as \code{x} with the density
#' (or log-density if \code{log=TRUE}).
#'
#' @examples
#' # Density at a point and on a vector
#' dllogis_ll(1, scale = 2, shape = 1.5)
#' dllogis_ll(c(0.5, 1, 3), scale = 2, shape = 1.5, log = TRUE)
#'
#' # Nonpositive x: 0 in density, -Inf in log-density
#' dllogis_ll(c(-1, 0), scale = 2, shape = 1.5)
#' dllogis_ll(0, scale = 2, shape = 1.5, log = TRUE)
#'
#' @seealso \code{\link{pllogis_ll}}
#' @family loglogistic
#' @export
dllogis_ll <- function(x, scale, shape, log = FALSE) {
  eta <- scale; xi <- shape
  if (!(eta > 0 && xi > 0)) return(rep(if (log) -Inf else 0, length(x)))
  x <- as.numeric(x)
  out <- rep(if (log) -Inf else 0, length(x))
  pos <- is.finite(x) & (x > 0)
  if (any(pos)) {
    z <- x[pos] / eta
    # log f = log(xi) - log(eta) + (xi - 1)*log(z) - 2*log1p(z^xi)
    logf <- log(xi) - log(eta) + (xi - 1) * log(z) - 2 * log1p(z^xi)
    out[pos] <- if (log) logf else exp(logf)
  }
  out
}
