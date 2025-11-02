#' Generalized Gamma density (Stacy form with f(t) = k * beta^(k*alpha) / Gamma(alpha) * t^(k*alpha-1) * exp(-(beta*t)^k))
#'
#' @description Density of the generalized gamma distribution with parameters
#' alpha (>0), beta (>0), kappa (>0):
#'   f(t) = [kappa * beta^(kappa*alpha) / Gamma(alpha)] * t^(kappa*alpha - 1) * exp(-(beta*t)^kappa),  t>0.
#'
#' @param x Numeric vector of quantiles (x > 0).
#' @param alpha Shape parameter (>0).
#' @param beta  Scale-rate parameter (>0) (appears as (beta*x)^kappa).
#' @param kappa Family/shape parameter (>0).
#' @param log Logical; if TRUE, return log-density.
#'
#' @return Numeric vector of (log-)densities of the same length as x.
#' @examples
#' dggamma(1, alpha=2, beta=1, kappa=0.5)
#' dggamma(c(0.5,1,2), alpha=3, beta=2, kappa=1.5, log=TRUE)
#'
#' @export
dggamma <- function(x, alpha, beta, kappa, log = FALSE) {
  x <- as.numeric(x)
  n <- length(x)

  # Parámetros fuera del dominio: devolver 0 (o -Inf en log) con la misma longitud
  if (!(alpha > 0 && beta > 0 && kappa > 0)) {
    return(rep(if (log) -Inf else 0, n))
  }

  # Sin datos: devolver vector vacío coherente
  if (n == 0) return(if (log) numeric(0) else numeric(0))

  # Plantilla de log-densidad: -Inf para x<=0 o no finitos
  logf <- rep(-Inf, n)

  # Índice vectorizado de valores válidos
  pos <- is.finite(x) & (x > 0)

  if (any(pos)) {
    xx <- x[pos]
    t  <- (beta * xx)^kappa
    # log f(t) = log kappa + (kappa*alpha) log beta - lgamma(alpha)
    #           + (kappa*alpha - 1) log t - (beta t)^kappa
    logf[pos] <- log(kappa) + (kappa * alpha) * log(beta) - lgamma(alpha) +
      (kappa * alpha - 1) * log(xx) - t
  }

  if (log) logf else exp(logf)
}
