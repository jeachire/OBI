#' Simple convergence diagnostics for random-walk Metropolis chains
#'
#' Computes basic, heuristic convergence diagnostics for a random-walk
#' Metropolis(-Hastings) Markov chain. The checks include:
#' \itemize{
#'   \item sanity checks on burn-in length and finite values;
#'   \item detection of numerical explosions via a maximum allowed absolute value;
#'   \item stability of the posterior mean, by comparing the first and second half
#'         of the post--burn-in segment;
#'   \item a reasonable range for the overall acceptance rate.
#' }
#' The goal is to flag clearly problematic chains in large simulation studies,
#' rather than to provide a formal convergence test.
#'
#' @param chain Numeric matrix with MCMC draws. Rows correspond to iterations,
#'   columns to parameters.
#' @param burnin Integer. Number of initial iterations to discard when computing
#'   diagnostics.
#' @param acc.rate Numeric scalar. Overall acceptance rate of the sampler,
#'   typically \code{accepted / n.iter}.
#' @param rel.tol Numeric scalar. Maximum allowed relative difference between
#'   the posterior means of the first and second half of the post--burn-in
#'   segment (default \code{0.05}).
#' @param max.abs Numeric scalar. Maximum allowed absolute value in the
#'   post--burn-in draws. Values above this threshold are interpreted as
#'   numerical explosions (default \code{1e6}).
#' @param acc.min Numeric scalar. Lower bound for a reasonable acceptance rate
#'   (default \code{0.05}).
#' @param acc.max Numeric scalar. Upper bound for a reasonable acceptance rate
#'   (default \code{0.7}).
#'
#' @return A list with components:
#' \describe{
#'   \item{converged}{Logical. \code{TRUE} if all checks are passed, \code{FALSE} otherwise.}
#'   \item{reason}{Character string describing the main reason for failure
#'     (e.g. \code{"burnin_too_large"}, \code{"non_finite_values"},
#'     \code{"exploding_values"}, \code{"means_or_accept_out_of_range"},
#'     or \code{"ok"} when \code{converged = TRUE}.}
#'   \item{rel.diff}{Numeric vector with the relative differences between the
#'     posterior means of the first and second half of the post--burn-in
#'     segment, for each parameter.}
#'   \item{acc.rate}{Numeric scalar. The acceptance rate supplied to the function.}
#' }
#'
#' @keywords internal
#' @noRd
diagnostic_RWmcmc <- function(chain,
                              burnin,
                              acc.rate,
                              rel.tol = 0.05,
                              max.abs = 1e6,
                              acc.min = 0.05,
                              acc.max = 0.7) {
  n.iter <- nrow(chain)

  # 1) Comprobaciones básicas
  if (burnin >= n.iter - 10L) {
    return(list(
      converged = FALSE,
      reason    = "burnin_too_large",
      rel.diff  = NA,
      acc.rate  = acc.rate
    ))
  }

  post <- chain[(burnin + 1L):n.iter, , drop = FALSE]

  if (any(!is.finite(post))) {
    return(list(
      converged = FALSE,
      reason    = "non_finite_values",
      rel.diff  = NA,
      acc.rate  = acc.rate
    ))
  }

  if (max(abs(post)) > max.abs) {
    return(list(
      converged = FALSE,
      reason    = "exploding_values",
      rel.diff  = NA,
      acc.rate  = acc.rate
    ))
  }

  # 2) Estabilidad de la media (primera mitad vs segunda mitad)
  mid <- floor(nrow(post) / 2)
  first_half  <- post[1:mid, , drop = FALSE]
  second_half <- post[(mid + 1):nrow(post), , drop = FALSE]

  mean1 <- colMeans(first_half)
  mean2 <- colMeans(second_half)

  rel.diff <- abs(mean2 - mean1) / (abs(mean2) + 1e-8)
  ok_means <- all(rel.diff < rel.tol)

  # 3) Rango razonable de aceptación
  ok_accept <- (acc.rate >= acc.min && acc.rate <= acc.max)

  converged <- ok_means && ok_accept

  list(
    converged = converged,
    reason    = if (converged) "ok" else "means_or_accept_out_of_range",
    rel.diff  = rel.diff,
    acc.rate  = acc.rate
  )
}
