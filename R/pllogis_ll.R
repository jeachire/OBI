#' @title Log-logistic cumulative distribution function
#' @name pllogis_ll
#' @description Computes the CDF of the log-logistic distribution for
#' \eqn{x \ge 0} under \code{scale = \eqn{\eta} > 0} and
#' \code{shape = \eqn{\xi} > 0}.
#'
#' @details
#' For \eqn{x>0}, the CDF is
#' \deqn{
#' F(x;\eta,\xi) \;=\; \frac{(x/\eta)^{\xi}}{1 + (x/\eta)^{\xi}}
#' \;=\; \left[1 + \left(\frac{x}{\eta}\right)^{-\xi}\right]^{-1}.
#' }
#' For \eqn{x \le 0}, \code{pllogis_ll} returns 0 in the lower tail
#' (\code{lower.tail = TRUE}) and 1 in the upper tail (\code{lower.tail = FALSE});
#' with \code{log.p = TRUE} it returns \code{-Inf} or \code{0}, respectively.
#'
#' @param x Numeric vector of quantiles (\eqn{x}).
#' @param scale Positive scale parameter \eqn{\eta > 0}.
#' @param shape Positive shape parameter \eqn{\xi > 0}.
#' @param lower.tail Logical; if \code{TRUE} (default) returns \eqn{P(X \le x)};
#' otherwise returns \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, returns probabilities on the log scale.
#'
#' @return A numeric vector of the same length as \code{x} with probabilities
#' (or log-probabilities if \code{log.p = TRUE}).
#'
#' @examples
#' # CDF and survival
#' pllogis_ll(1, scale = 2, shape = 1.5)                   # F(1)
#' pllogis_ll(1, scale = 2, shape = 1.5, lower.tail = FALSE) # S(1)
#'
#' # Log probabilities (useful in censored log-likelihoods)
#' pllogis_ll(c(0.5, 1, 3), scale = 2, shape = 1.5,
#'            lower.tail = FALSE, log.p = TRUE)
#'
#' # Behavior at x <= 0
#' pllogis_ll(c(-1, 0), scale = 2, shape = 1.5)                 # 0
#' pllogis_ll(0, scale = 2, shape = 1.5, lower.tail = FALSE)    # 1
#' pllogis_ll(0, scale = 2, shape = 1.5, lower.tail = FALSE, log.p = TRUE)  # 0
#'
#' @seealso \code{\link{dllogis_ll}}
#' @family loglogistic
#' @export
pllogis_ll <- function(x, scale, shape, lower.tail = TRUE, log.p = FALSE) {
  eta <- scale; xi <- shape
  if (!(eta > 0 && xi > 0)) {
    return(rep(if (log.p) -Inf else if (lower.tail) 0 else 1, length(x)))
  }
  x <- as.numeric(x)
  out <- rep(if (log.p) -Inf else if (lower.tail) 0 else 1, length(x))
  pos <- is.finite(x) & (x > 0)
  if (any(pos)) {
    z <- x[pos] / eta
    a <- z^xi
    if (lower.tail) {
      # F = a / (1 + a); log F = log(a) - log1p(a)
      val <- if (log.p) (log(a) - log1p(a)) else (a / (1 + a))
    } else {
      # S = 1/(1 + a); log S = -log1p(a)
      val <- if (log.p) (-log1p(a)) else (1 / (1 + a))
    }
    out[pos] <- val
  }
  # Para x<=0: F=0 (cola inferior) y S=1 (cola superior) ya está manejado por defaults
  out
}
