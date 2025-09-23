#' Bayesian Fit of Parametric Lifetime Distributions
#'
#' Fits a parametric distribution to (possibly censored) data under an objective prior,
#' and samples from the posterior distribution using MCMC. The maximum likelihood
#' estimator (MLE) is computed and used as the starting value for the MCMC sampler.
#'
#' @param x Numeric vector. Observed data (event times or counts).
#' @param dist Character string. Distribution family to fit.
#'   Supported options: \code{"normal"}, \code{"poisson"}, \code{"gamma"},
#'   \code{"exponential"}, \code{"weibull"}.
#' @param delta Integer or logical vector of the same length as \code{x}.
#'   Censoring indicator: \code{1} for observed values, \code{0} for right-censored.
#'   If \code{NULL}, all observations are assumed to be uncensored.
#' @param prior Character string or function. Prior distribution to use.
#'   Built-in options: \code{"Jeffreys"}, \code{"Reference"}, \code{"MDI"}, \code{"Uniform"}.
#'   Alternatively, the user may provide a custom prior function.
#' @param method Character string. MCMC method to use.
#'   \code{"MCMC"} runs a random walk Metropolis-Hastings;
#'   \code{"MCMC_I"} runs Metropolis-Hastings with imputation of right-censored data.
#' @param n.iter Integer. Number of total MCMC iterations (default 5000).
#' @param burnin Integer. Number of burn-in iterations (default 1000).
#' @param thin Integer. Thinning interval, i.e., keep every \code{thin}-th sample (default 1).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{chain}{Matrix of sampled parameter values after burn-in and thinning.}
#'   \item{summary}{Data frame with posterior mean, SD, median, and quantiles (2.5\%, 25\%, 75\%, 97.5\%).}
#'   \item{acc.rate}{Acceptance rate of the sampler.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The function automatically determines which parameters are restricted to the
#' positive domain, applies a log-transformation where needed, and adjusts the
#' posterior with the Jacobian. For censored data, the \code{"MCMC_I"} method
#' imputes censored values at each iteration.
#'
#' @examples
#' ## Example: Fit a Weibull distribution with right censoring
#' set.seed(123)
#' n <- 40
#' t.true <- rweibull(n, shape=2, scale=1)
#' censor.limit <- 1.5
#' x <- pmin(t.true, censor.limit)
#' delta <- as.integer(t.true <= censor.limit)
#'
#' # Fit model
#' fit <- fitBayes(
#'   x = x, dist = "weibull", delta = delta,
#'   prior = "Jeffreys", method = "MCMC_I",
#'   n.iter = 2000, burnin = 500
#' )
#'
#' fit$summary
#'
#' @export
fitBayes <- function(x, dist, delta = NULL, prior = "Jeffreys",
                     method = "MCMC", n.iter = 5000, burnin=1000,
                     thin=1) {

  # 1. Validaciones
  stopifnot(is.numeric(x))
  stopifnot(dist %in% c("normal","poisson","gamma","exponential","weibull"))

  if(is.null(delta)) delta <- rep(1, length(x))  # todos observados si no se indica censura
  stopifnot(length(delta) == length(x))

  # 2. Log-verosimilitud (adaptada para censura)
  loglik <- getLogLikelihood(dist, x, delta)
  positive <- getPositive(dist)

  # 3. Prior objetiva
  priorFunc <- getPrior(dist, prior)

  # 4. Posterior
  logpost <- function(theta) loglik(theta) + priorFunc(theta)

  logPost <- function(theta, t) {
    logLik <- getLogLikelihood(dist, t)   # logLik ahora es una función de theta
    logLik(theta) + priorFunc(theta)      # evaluamos theta
  }



  # 5. Theta inicial
  p <- length(positive)
  theta0 <- ifelse(positive, 1, 0)
  lower <- ifelse(positive, 1e-6, -Inf)
  upper <- rep(Inf,p)
  theta.init <- optim(par = theta0,
               fn = function(theta) -loglik(theta),
               method = "L-BFGS-B",
               lower = lower,
               upper = upper)$par

  # 6. Muestra posterior
  if(method == "MCMC"){
    postSample <- logRWmcmc(logPost = logpost,
                            theta.init = theta.init,
                            positive = positive,
                            Sigma = NULL,
                            n.iter = n.iter,
                            burnin = burnin,
                            thin = thin)

    } else if (method == "MCMC_I"){
      postSample <- censoredRWmcmc(x = x,
                              delta = delta,
                              dist = dist,
                              logPost = logPost,
                              theta.init = theta.init,
                              positive = positive,
                              Sigma = NULL,
                              n.iter = n.iter,
                              burnin = burnin,
                              thin = thin)
    } else {
      stop("Solo 'MCMC' y 'MCMC_I' está implementado actualmente.")
    }

  # 7. Estimaciones resumen
  est <- postSample$summary
  posterior <- postSample$chain
  rate <- postSample$acc.rate

  # 8. Salida
  return(list(
    summary = est,
    posterior = posterior,
    rate = rate,
    call = match.call()
  ))
}

