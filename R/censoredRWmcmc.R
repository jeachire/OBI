#' Random Walk Metropolis-Hastings with Right-Censoring Imputation
#'
#' Runs a Metropolis-Hastings MCMC algorithm with Gaussian random-walk proposals,
#' imputing right-censored observations at each iteration. Parameters constrained
#' to be positive are handled via a log-transformation.
#'
#' @param x Numeric vector. Observed times (may include censored observations).
#' @param delta Integer or logical vector of the same length as \code{x}.
#'   Indicator of censoring: \code{1} for observed events, \code{0} for right-censored.
#' @param dist Character string. Distribution used for imputing censored values.
#'   Currently supports \code{"weibull"} and \code{"exponential"}.
#' @param logPost Function. Log-posterior function to evaluate.
#'   Must take arguments \code{theta} (numeric vector of parameters)
#'   and \code{t.imp} (data vector with imputed values).
#' @param theta.init Numeric vector. Initial parameter values.
#' @param positive Logical vector of the same length as \code{theta.init}.
#'   Indicates which parameters are constrained to be positive. Default: all \code{FALSE}.
#' @param Sigma Covariance matrix for the Gaussian proposal distribution.
#'   If \code{NULL}, a diagonal matrix with variance 0.1 is used.
#' @param n.iter Integer. Total number of MCMC iterations (default 5000).
#' @param burnin Integer. Number of initial iterations to discard (default 1000).
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
#' At each iteration, censored observations are imputed from the conditional
#' distribution of the failure time given censoring. The sampler then performs
#' a Metropolis-Hastings step using the imputed dataset. The Jacobian adjustment
#' for positive parameters is automatically included.
#'
#' @examples
#' # Simulated censored data (Weibull; parameterization = (rate, shape))
#' set.seed(123)
#' n <- 30
#' shape <- 2; rate <- 1
#' t.true <- rweibull(n, shape = shape, scale = 1/rate)
#' censor.limit <- 1.2
#' x <- pmin(t.true, censor.limit)
#' delta <- as.integer(t.true <= censor.limit)
#'
#' # Log-posterior (up to a constant) in (rate, shape)
#' logPost <- function(theta, t.imp){
#'   rate  <- theta[1]; shape <- theta[2]
#'   scale <- 1 / rate
#'   sum(dweibull(t.imp, shape = shape, scale = scale, log = TRUE))
#' }
#'
#' out <- censoredRWmcmc(
#'   x = x, delta = delta, dist = "weibull",
#'   logPost = logPost, theta.init = c(rate = 1, shape = 2),
#'   positive = c(TRUE, TRUE), n.iter = 2000
#' )
#' str(out$summary)
#'
#' @importFrom MASS mvrnorm
#' @export
censoredRWmcmc <- function(x, delta, dist, logPost, theta.init, positive,
                           Sigma = NULL, n.iter=5000, burnin=1000,
                           thin=1){

  n <- length(x)
  p <- length(theta.init)

  if(is.null(Sigma)) Sigma <- 0.1 * diag(p)
  if(is.null(positive)) positive <- rep(FALSE, p)

  # Transformaciones log para parámetros positivos
  transform_fwd <- function(theta){
    phi <- theta
    phi[positive] <- log(theta[positive])
    phi
  }
  transform_inv <- function(phi){
    theta <- phi
    theta[positive] <- exp(phi[positive])
    theta
  }

  # Jacobiano logarítmico
  logJacobian <- function(theta){
    sum(log(theta[positive]))
  }

  # Inicialización
  theta.current <- theta.init
  phi.current <- transform_fwd(theta.current)
  chain <- matrix(NA, nrow=n.iter, ncol=p)
  acc <- 0

  for(j in 1:n.iter){
    # --- 1. Imputar censura ---
    t.imp <- imputeRightCensor(x=x, delta=delta, theta=theta.current, dist=dist)

    # --- 2. Propuesta theta* ---
    phi.prop <- mvrnorm(1, mu=phi.current, Sigma=Sigma)
    theta.prop <- transform_inv(phi.prop)

    # --- 3. Ratio de Metropolis ---
    logpost.prop <- logPost(theta.prop, t.imp) + logJacobian(theta.prop)
    logpost.curr <- logPost(theta.current, t.imp) + logJacobian(theta.current)

    log.alpha <- logpost.prop - logpost.curr

    if(log(runif(1)) < log.alpha){
      theta.current <- theta.prop
      phi.current <- phi.prop
      acc <- acc + 1
    }

    chain[j,] <- theta.current
  }

  # --- Postprocesamiento ---
  keep <- seq(from=burnin+1, to=n.iter, by=thin)
  chain.kept <- chain[keep,,drop=FALSE]

  summary <- data.frame(
    Mean   = apply(chain.kept, 2, mean),
    SD     = apply(chain.kept, 2, sd),
    Median = apply(chain.kept, 2, median),
    Q2.5   = apply(chain.kept, 2, quantile, 0.025),
    Q25    = apply(chain.kept, 2, quantile, 0.25),
    Q75    = apply(chain.kept, 2, quantile, 0.75),
    Q97.5  = apply(chain.kept, 2, quantile, 0.975)
  )

  return(list(
    chain = chain.kept,
    summary = summary,
    acc.rate = acc/n.iter,
    call = match.call()
  ))
}
