#' @title Bayesian Fit of Parametric Lifetime Distributions (with optional right-censoring)
#'
#' @description
#' Fits parametric lifetime models under an objective (or user-chosen) prior by
#' constructing the log-likelihood and log-prior, then running a Metropolis–Hastings
#' MCMC sampler. When right-censored data are present, the function can perform
#' data augmentation by imputing censored observations at each MCMC iteration.
#' Returns posterior draws, posterior summaries, and a frequentist MLE used for initialization.
#'
#' @param x Numeric vector of data (event times or counts). May include censored observations.
#' @param dist Character string naming the distribution. Supported options:
#' \itemize{
#'   \item \code{"normal"}: \code{theta = c(mean, variance)}
#'   \item \code{"poisson"}: \code{theta = lambda}
#'   \item \code{"exponential"}: \code{theta = c(rate)}
#'   \item \code{"weibull"}: \code{theta = c(rate, shape)} \emph{(note: base R uses \code{(shape, scale)};
#'         internally \code{scale = 1/rate})}
#'   \item \code{"gamma"}: \code{theta = c(shape, rate)}
#'   \item \code{"generalized gamma"} (Stacy): \code{theta = c(alpha, beta, kappa)}
#'   \item \code{"lognormal"}: \code{theta = c(meanlog, varlog)} \emph{(so \code{sdlog = sqrt(varlog)})}
#'   \item \code{"loglogistic"}: \code{theta = c(scale, shape)}
#' }
#' @param delta Optional numeric/logical vector the same length as \code{x} with censoring indicators:
#' \code{1} for observed events, \code{0} for right-censored. If \code{NULL}, all observations are treated as observed.
#' @param prior Prior choice for \code{theta}. Either a character string naming a built-in objective prior
#' (e.g., \code{"Jeffreys"}, \code{"reference"}, \code{"MDI"})—handled internally by \code{getPrior()}—or a user-supplied
#' function returning a log-prior in \code{theta}-space.
#' @param method Character string. Sampling method:
#' \itemize{
#'   \item \code{"MCMC"} – standard random-walk MH on the transformed parameter space.
#'   \item \code{"MCMC_I"} – MH with right-censoring imputation (data augmentation) at each iteration.
#' }
#' @param n.iter Integer. Total number of MCMC iterations (default \code{5000}).
#' @param burnin Integer. Number of initial iterations to discard (default \code{1000}).
#' @param thin Integer. Keep every \code{thin}-th draw (default \code{1}).
#' @param verbose Logical; if \code{TRUE}, prints basic initialization info (default \code{TRUE}).
#'
#' @details
#' \strong{Likelihood and prior:}
#' The function builds a log-likelihood via \code{getLogLikelihood(dist, x, delta)} and a log-prior via
#' \code{getPrior(dist, prior)} under the parameterizations listed in \emph{dist}.
#'
#' \strong{Transformations and Jacobian:}
#' Parameters restricted to \((0,\infty)\) are automatically mapped to the real line through a log-transform.
#' The sampler accounts for the change of variables by adding the log-Jacobian
#' \eqn{\sum_{j:\,\theta_j>0}\log(\theta_j)} to the target density.
#'
#' \strong{Censoring and data augmentation (\code{method = "MCMC_I"}):}
#' For right-censored observations (\code{delta = 0}), the algorithm imputes at each iteration from
#' the conditional distribution \eqn{T \mid T > c} using inverse-CDF sampling. The updated complete data are then
#' used in the MH step. Implementations rely only on \pkg{stats} functions (and simple helpers for log-logistic
#' and generalized gamma).
#'
#' \strong{Initialization:}
#' The sampler is initialized at a frequentist MLE obtained by optimizing the log-likelihood on the transformed scale.
#' The inverse observed Hessian at the MLE is used as the default proposal covariance on the transformed scale, with
#' fallbacks if the Hessian is not numerically invertible.
#'
#' \strong{Parameterizations (summary):}
#' \tabular{ll}{
#'   Normal \tab \code{(mean, variance)} \cr
#'   Poisson \tab \code{lambda} \cr
#'   Exponential \tab \code{(rate)} \cr
#'   Weibull \tab \code{(rate, shape)} \; (\code{scale = 1/rate} for base R calls) \cr
#'   Gamma \tab \code{(shape, rate)} \cr
#'   Generalized gamma (Stacy) \tab \code{(alpha, beta, kappa)} with \eqn{Y = (\beta T)^\kappa \sim \mathrm{Gamma}(\alpha, 1)} \cr
#'   Lognormal \tab \code{(meanlog, varlog)} \cr
#'   Log-logistic \tab \code{(scale, shape)}
#' }
#'
#' \strong{Priors:}
#' Built-in objective priors (Jeffreys, reference, MDI) follow the package’s documented conventions and may be
#' improper; posterior propriety depends on the data and model. The user may supply a custom log-prior function
#' if desired.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{chain}}{Matrix of posterior draws for \code{theta} (after burn-in and thinning).}
#'   \item{\code{summary}}{Data frame with posterior mean, SD, median, and central quantiles (2.5\%, 25\%, 75\%, 97.5\%).}
#'   \item{\code{acc.rate}}{Overall Metropolis–Hastings acceptance rate.}
#'   \item{\code{diagnostics}}{List with simple convergence diagnostics returned by the underlying
#'         Metropolis sampler (see \code{\link{logRWmcmc}} and \code{\link{censoredRWmcmc}}).}
#'   \item{\code{theta_mle}}{Frequentist MLE used for initialization (on the natural \code{theta}-scale).}
#'   \item{\code{Sigma_init_phi}}{Initial proposal covariance on the transformed scale (inverse Hessian at the MLE when available).}
#'   \item{\code{call}}{Matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Example: Weibull with right-censoring, objective prior, MH with imputation
#' set.seed(123)
#' n <- 50
#' shape <- 2
#' rate  <- 1                   # scale = 1 / rate = 1
#' t_true <- rweibull(n, shape = shape, scale = 1 / rate)
#' c_lim  <- 1.2
#' x      <- pmin(t_true, c_lim)
#' delta  <- as.integer(t_true <= c_lim)  # 1 = observed, 0 = censored
#'
#' fit <- fitBayes(
#'   x = x,
#'   dist = "weibull",                   # parameterization: theta = c(rate, shape)
#'   delta = delta,
#'   prior = "Jeffreys",
#'   method = "MCMC_I",                  # impute censored data each iteration
#'   n.iter = 4000, burnin = 1000, thin = 2
#' )
#'
#' fit$summary
#' head(fit$chain)
#' fit$acc.rate
#' fit$theta_mle
#' }
#'
#' @seealso \code{\link{getLogLikelihood}}, \code{\link{getPrior}}, \code{\link{imputeRightCensor}},
#' \code{\link{logRWmcmc}}, \code{\link{censoredRWmcmc}}
#'
#' @note
#' Be mindful of parameterizations. In particular, \code{weibull} uses \code{(rate, shape)} in this package
#' and converts to base R’s \code{(shape, scale)} internally via \code{scale = 1/rate}.
#' Some objective priors are improper; ensure posterior propriety for your dataset/model.
#'
#' @references
#' Bernardo, J. M., & Smith, A. F. M. (2000). \emph{Bayesian Theory}. Wiley.
#' Berger, J. O., Bernardo, J. M., & Sun, D. (2009). The formal definition of reference priors.
#' \emph{Annals of Statistics}, 37(2), 905–938.
#' Kleiber, C., & Kotz, S. (2003). \emph{Statistical Size Distributions in Economics and Actuarial Sciences}. Wiley.
#'
#' @export
fitBayes <- function(x, dist, delta = NULL, prior = "Jeffreys",
                     method = "MCMC", n.iter = 5000, burnin = 1000, thin = 1,
                     verbose = TRUE) {

  ## ---------- 0) Basic checks ----------
  if (!is.numeric(x)) stop("'x' must be numeric.")
  if (is.null(delta)) delta <- rep(1L, length(x))
  if (length(delta) != length(x)) stop("'delta' must have the same length as 'x'.")

  supported <- c("normal","poisson","gamma","exponential","weibull",
                 "loglogistic","lognormal","generalized gamma")
  if (!dist %in% supported) stop("Unsupported distribution: ", dist)

  ## ---------- 1) Build log-likelihood and prior (on theta / natural scale) ----------
  loglik_fun <- getLogLikelihood(dist, x, delta)     # returns function(theta)
  positive   <- getPositive(dist)                    # logical vector by your conventions
  p <- length(positive)

  prior_fun_user <- getPrior(dist, prior)            # returns log-prior function(theta)
  if (!is.function(prior_fun_user)) stop("getPrior must return a function.")

  # If you keep the heuristic, leave as is; otherwise assume it's a log-prior already:
  prior_fun_log <- function(theta) prior_fun_user(theta)

  ## ---------- 2) Transformations phi <-> theta ----------
  idx_pos <- which(positive)
  idx_un  <- which(!positive)

  to_phi <- function(theta) {
    phi <- numeric(p)
    if (length(idx_pos)) phi[idx_pos] <- log(theta[idx_pos])
    if (length(idx_un))  phi[idx_un]  <- theta[idx_un]
    phi
  }
  from_phi <- function(phi) {
    theta <- numeric(p)
    if (length(idx_pos)) theta[idx_pos] <- exp(phi[idx_pos])
    if (length(idx_un))  theta[idx_un]  <- phi[idx_un]
    theta
  }
  logJ_phi <- function(phi) {
    if (!length(idx_pos)) return(0)
    sum(phi[idx_pos])  # since theta_pos = exp(phi_pos)
  }

  ## ---------- 3) MLE (in phi-space) for initialization ----------
  opt <- tryCatch(
    optim(par = rep(0, p),
          fn = function(phi) {
            th <- from_phi(phi)
            val <- loglik_fun(th)
            if (!is.finite(val)) return(1e12)  # crude guard
            -val
          },
          method = "BFGS", hessian = TRUE, control = list(maxit = 1000)),
    error = function(e) NULL
  )


  # --- Robust initialization in phi-space (MLE/MAP) and Sigma_phi0 ---
  if (!is.null(opt) && is.numeric(opt$par) && isTRUE(opt$convergence == 0)) {
    # Use optimizer result
    phi_mle   <- opt$par
    theta_mle <- from_phi(phi_mle)

    # Try to form proposal covariance from inverse Hessian; fallback to diag(0.2, p)
    Sigma_phi0 <- tryCatch({
      S <- solve(opt$hessian)
      # Symmetrize defensively
      S <- (S + t(S)) / 2
      # If any non-finite elements, force fallback
      if (any(!is.finite(S))) stop("non-finite covariance")
      S
    }, error = function(e) diag(0.2, p))

    # Final guard: if still non-finite, fallback
    if (any(!is.finite(Sigma_phi0))) {
      Sigma_phi0 <- diag(0.2, p)
    }

  } else {
    # No convergence or invalid opt -> neutral initialization
    phi_mle   <- rep(0, p)        # from_phi will map positives -> 1, unbounded -> 0
    theta_mle <- from_phi(phi_mle)
    Sigma_phi0 <- diag(0.2, p)
    warning("optim did not converge; using neutral initial values (phi=0, Sigma=0.2*I).")
  }

  if (verbose) message("Initial theta (natural): ",
                       paste0(round(theta_mle, 4), collapse = ", "))


  ## ---------- 4) Define targets for samplers ----------
  # (a) Target on THETA: loglik + logprior (NO Jacobian here)
  logPost_theta <- function(theta) {
    ll <- loglik_fun(theta)
    lp <- prior_fun_log(theta)
    if (!is.finite(ll) || !is.finite(lp)) return(-Inf)
    ll + lp
  }

  # (b) log-posterior en theta usando los datos imputados t.imp

  logPost_theta_with_timp <- function(theta, t.imp) {
    ll_fun <- getLogLikelihood(
      dist  = dist,
      x     = t.imp,
      delta = rep(1L, length(t.imp))
    )
    ll <- ll_fun(theta)
    lp <- prior_fun_log(theta)
    if (!is.finite(ll) || !is.finite(lp)) return(-Inf)
    ll + lp
  }

  ## ---------- 5) Run sampler ----------
  if (method == "MCMC") {
    res <- logRWmcmc(
      logPost     = logPost_theta,
      theta.init  = theta_mle,
      positive    = positive,
      Sigma       = Sigma_phi0,
      n.iter      = n.iter, burnin = burnin, thin = thin
    )
    chain_theta <- res$chain
    acc_rate    <- res$acc.rate
    diag_obj    <- if (!is.null(res$diagnostics)) res$diagnostics else NULL

  } else if (method == "MCMC_I") {

    res <- censoredRWmcmc(
      x          = x,
      delta      = delta,
      dist       = dist,
      logPost    = logPost_theta_with_timp,
      theta.init = theta_mle,
      positive   = positive,
      Sigma      = Sigma_phi0,
      n.iter     = n.iter, burnin = burnin, thin = thin
    )
    chain_theta <- res$chain
    acc_rate    <- res$acc.rate
    diag_obj    <- if (!is.null(res$diagnostics)) res$diagnostics else NULL

  } else {
    stop("Only 'MCMC' and 'MCMC_I' are implemented.")
  }

  ## ---------- 6) Summaries ----------

  if (verbose && !is.null(diag_obj) && !isTRUE(diag_obj$converged)) {
    warning("MCMC diagnostics indicate potential non-convergence (reason = '",
            diag_obj$reason, "').", call. = FALSE)
  }

  summary_df <- data.frame(
    Mean   = apply(chain_theta, 2, mean),
    SD     = apply(chain_theta, 2, sd),
    Median = apply(chain_theta, 2, median),
    Q2.5   = apply(chain_theta, 2, quantile, 0.025),
    Q25    = apply(chain_theta, 2, quantile, 0.25),
    Q75    = apply(chain_theta, 2, quantile, 0.75),
    Q97.5  = apply(chain_theta, 2, quantile, 0.975),
    check.names = FALSE
  )

  list(
    chain          = chain_theta,
    summary        = summary_df,
    acc.rate       = acc_rate,
    diagnostics    = diag_obj,
    call           = match.call(),
    theta_mle      = theta_mle,
    Sigma_init_phi = Sigma_phi0
  )
}
