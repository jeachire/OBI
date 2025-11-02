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
  # Coerción y longitud
  x <- as.numeric(x)
  n <- length(x)

  # Chequeo de dominio de parámetros: devolvemos NA del mismo largo
  if (!(alpha > 0 && beta > 0 && kappa > 0)) {
    return(rep(if (log.p) NaN else NA_real_, n))
  }

  # Vector resultado (CDF o log-CDF, según log.p)
  res <- rep(NA_real_, n)

  # Máscaras útiles
  fin  <- is.finite(x)
  neg0 <- fin & (x <= 0)        # incluye x == 0 y x < 0
  pos  <- fin & (x > 0)
  pinf <- is.infinite(x) & (x > 0)  # +Inf
  ninf <- is.infinite(x) & (x < 0)  # -Inf (trátalo como x<=0)

  # x <= 0  --> F(x)=0   ; S(x)=1
  if (any(neg0 | ninf)) {
    idx <- (neg0 | ninf)
    if (log.p) {
      res[idx] <- if (lower.tail) -Inf else 0     # log(0)=-Inf, log(1)=0
    } else {
      res[idx] <- if (lower.tail) 0    else 1
    }
  }

  # x == +Inf  --> F(x)=1 ; S(x)=0
  if (any(pinf)) {
    if (log.p) {
      res[pinf] <- if (lower.tail) 0   else -Inf  # log(1)=0, log(0)=-Inf
    } else {
      res[pinf] <- if (lower.tail) 1   else 0
    }
  }

  # x > 0 y finito: usar transformación (beta*x)^kappa ~ Gamma(shape=alpha, rate=1)
  if (any(pos)) {
    t <- (beta * x[pos])^kappa
    res[pos] <- stats::pgamma(t, shape = alpha, rate = 1,
                              lower.tail = lower.tail, log.p = log.p)
  }

  # x = NA/NaN quedan como NA (propagan NA)
  res
}

