#' Simple convergence diagnostics for random-walk Metropolis chains
#'
#' Computes basic convergence diagnostics for a random-walk
#' Metropolis(-Hastings) Markov chain. The checks include:
#' \itemize{
#'   \item sanity checks on burn-in length and finite values;
#'   \item detection of numerical explosions via a maximum allowed absolute value;
#'   \item a reasonable range for the overall acceptance rate;
#'   \item a Geweke-type diagnostic for each parameter based on the
#'         \pkg{coda} implementation.
#' }
#' The goal is to flag clearly problematic chains in large simulation studies,
#' rather than to provide a fully exhaustive convergence analysis.
#'
#' @param chain Numeric matrix with MCMC draws. Rows correspond to iterations,
#'   columns to parameters.
#' @param burnin Integer. Number of initial iterations to discard when computing
#'   diagnostics.
#' @param acc.rate Numeric scalar. Overall acceptance rate of the sampler,
#'   typically \code{accepted / n.iter}.
#' @param max.abs Numeric scalar. Maximum allowed absolute value in the
#'   post--burn-in draws. Values above this threshold are interpreted as
#'   numerical explosions (default \code{1e6}).
#' @param acc.min Numeric scalar. Lower bound for a reasonable acceptance rate
#'   (default \code{0.05}).
#' @param acc.max Numeric scalar. Upper bound for a reasonable acceptance rate
#'   (default \code{0.7}).
#' @param z.threshold Numeric scalar. Threshold for the absolute value of the
#'   Geweke \eqn{Z}-scores. Values below this threshold are interpreted as
#'   indicating no strong evidence against convergence (default \code{2}).
#'
#' @return A list with components:
#' \describe{
#'   \item{converged}{Logical. \code{TRUE} if all checks are passed, \code{FALSE} otherwise.}
#'   \item{reason}{Character string describing the main reason for failure
#'     (e.g. \code{"burnin_too_large"}, \code{"non_finite_values"},
#'     \code{"exploding_values"}, \code{"acceptance_out_of_range"},
#'     \code{"geweke_fail"}, or \code{"ok"} when \code{converged = TRUE}.}
#'   \item{geweke.z}{Numeric vector with the Geweke \eqn{Z}-scores for each parameter.}
#'   \item{acc.rate}{Numeric scalar. The acceptance rate supplied to the function.}
#' }
#'
#' @keywords internal
#' @noRd
#' @importFrom coda mcmc geweke.diag
diagnostic_RWmcmc <- function(chain,
                              burnin,
                              acc.rate,
                              max.abs    = 1e6,
                              acc.min    = 0.05,
                              acc.max    = 0.7,
                              z.threshold = 2) {

  n.iter <- nrow(chain)

  ## 1) Comprobaciones básicas sobre burn-in y longitud efectiva
  if (burnin >= n.iter - 10L) {
    return(list(
      converged = FALSE,
      reason    = "burnin_too_large",
      geweke.z  = NA_real_,
      acc.rate  = acc.rate
    ))
  }

  post <- chain[(burnin + 1L):n.iter, , drop = FALSE]

  ## 2) Valores finitos
  if (any(!is.finite(post))) {
    return(list(
      converged = FALSE,
      reason    = "non_finite_values",
      geweke.z  = NA_real_,
      acc.rate  = acc.rate
    ))
  }

  ## 3) Explosiones numéricas
  if (max(abs(post)) > max.abs) {
    return(list(
      converged = FALSE,
      reason    = "exploding_values",
      geweke.z  = NA_real_,
      acc.rate  = acc.rate
    ))
  }

  ## 4) Rango razonable de aceptación global
  ok_accept <- (acc.rate >= acc.min && acc.rate <= acc.max)
  if (!ok_accept) {
    return(list(
      converged = FALSE,
      reason    = "acceptance_out_of_range",
      geweke.z  = NA_real_,
      acc.rate  = acc.rate
    ))
  }

  ## 5) Diagnóstico clásico de Geweke (una cadena)
  #    Trabajamos parámetro a parámetro para evitar depender de la estructura
  #    interna multivariada del objeto 'geweke.diag'.
  p <- ncol(post)
  z_vals <- numeric(p)

  for (j in seq_len(p)) {
    mcj  <- coda::mcmc(post[, j])
    gj   <- coda::geweke.diag(mcj)
    z_j  <- as.numeric(gj$z)
    z_vals[j] <- z_j
  }

  # Comprobamos que no haya NAs y que todos los |Z| sean menores al umbral
  if (any(!is.finite(z_vals))) {
    return(list(
      converged = FALSE,
      reason    = "geweke_fail",
      geweke.z  = z_vals,
      acc.rate  = acc.rate
    ))
  }

  ok_geweke <- all(abs(z_vals) < z.threshold)

  converged <- ok_geweke

  list(
    converged = converged,
    reason    = if (converged) "ok" else "geweke_fail",
    geweke.z  = z_vals,
    acc.rate  = acc.rate
  )
}
