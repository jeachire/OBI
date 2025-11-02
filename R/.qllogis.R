.qllogis <- function(p, scale, shape) {
  eta <- scale; xi <- shape
  if (!(eta > 0 && xi > 0)) stop("Log-logistic requires scale > 0 and shape > 0.")
  # guard p in (0,1)
  p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
  eta * (p / (1 - p))^(1 / xi)
}
