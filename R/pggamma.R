#' Generalized Gamma CDF (Stacy form)
#'
#' @description Cumulative distribution function of the generalized gamma with
#' parameters alpha (>0), beta (>0), kappa (>0):
#'   F(t) = P(Y <= (beta*t)^kappa),  where Y ~ Gamma(shape=alpha, rate=1).
#'
#' @param x Numeric vector of quantiles (x >= 0).
#' @param alpha Shape parameter (>0).
#' @param beta  Scale-rate parameter (>0).
#' @param kappa Family/shape parameter (>0).
#' @param lower.tail Logical; if TRUE (default), returns P(X <= x), else P(X > x).
#' @param log.p Logical; if TRUE, returns log probabilities.
#'
#' @return Numeric vector of (log-)probabilities of the same length as x.
#' @examples
#' pggamma(1, alpha=2, beta=1, kappa=0.5)
#' pggamma(c(0, 0.5, 1, 2), alpha=3, beta=2, kappa=1.5, lower.tail=FALSE, log.p=TRUE)
#'
#' @export
pggamma <- function(x, alpha, beta, kappa, lower.tail = TRUE, log.p = FALSE) {
  if (!(alpha > 0 && beta > 0 && kappa > 0)) {
    return(rep(if (log.p) -Inf else if (lower.tail) 0 else 1, length(x)))
  }
  x <- as.numeric(x)
  # t = (beta*x)^kappa; usa pmax(x,0) para robustez en x<=0
  t <- (pmax(x, 0) * beta) ^ kappa
  stats::pgamma(t, shape = alpha, rate = 1, lower.tail = lower.tail, log.p = log.p)
}
