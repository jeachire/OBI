.qquantile_ggamma <- function(p, alpha, beta, kappa) {
  if (!(alpha > 0 && beta > 0 && kappa > 0))
    stop("Generalized gamma requires alpha > 0, beta > 0, kappa > 0.")
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  y <- qgamma(p, shape = alpha, rate = 1)
  (y)^(1 / kappa) / beta
}
