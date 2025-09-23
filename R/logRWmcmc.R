#' Random Walk Metropolis-Hastings with Log-Transformation
#'
#' Runs a Metropolis-Hastings MCMC algorithm with Gaussian random-walk proposals.
#' Parameters constrained to be positive are handled via a log-transformation.
#' The function also returns posterior summaries.
#'
#' @param logPost Function. The log-posterior function to evaluate.
#'   Must take a numeric vector \code{theta} as input and return a scalar.
#' @param theta.init Numeric vector. Initial values for the parameters.
#' @param positive Logical vector of the same length as \code{theta.init}.
#'   Indicates which parameters are constrained to be positive.
#' @param Sigma Covariance matrix for the Gaussian proposal distribution.
#'   If \code{NULL}, a diagonal matrix with variance 0.1 is used.
#' @param n.iter Integer. Total number of MCMC iterations (default 5000).
#' @param burnin Integer. Number of initial iterations to discard (default 1000).
#' @param thin Integer. Thinning interval, i.e., keep every \code{thin}-th sample (default 1).
#'
#' @return A list with the following elements:
#' \describe{
#' \item{chain}{Matrix of sampled parameter values after burn-in and thinning.}
#' \item{summary}{Data frame with posterior mean, SD, median, and quantiles (2.5\%, 25\%, 75\%, 97.5\%).}
#' \item{acc.rate}{Acceptance rate of the sampler.}
#' \item{call}{The matched call.}
#' }
#'
#' @details
#' The algorithm uses a Gaussian random-walk proposal in the transformed
#' parameter space. Parameters constrained to be positive are mapped via the
#' logarithm to the real line. The Jacobian adjustment is automatically added
#' to the log-posterior to ensure correct sampling.
#'
#' @examples
#' # Example: Normal likelihood with known variance
#' logPost <- function(theta){
#'   mu <- theta[1]
#'   sum(dnorm(rnorm(50, 2, 1), mean = mu, sd = 1, log = TRUE))
#' }
#' out <- logRWmcmc(logPost, theta.init = c(0), positive = c(FALSE), n.iter=2000)
#' str(out$summary)
#'
#' @importFrom MASS mvrnorm
#' @export
logRWmcmc <- function(logPost, theta.init, positive,
                      Sigma=NULL, n.iter=5000, burnin=1000,
                    thin=1){

  p <- length(theta.init)

  # --- Construir matriz de propuesta ---
  if(is.null(Sigma)) {
    Sigma <- (0.1) * diag(p)
  }

  # --- Transformaciones ---
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
  logJacobian <- function(theta){
    if(any(theta[positive] <= 0)) return(-Inf)
    sum(log(theta[positive]))
  }

  # --- Inicialización ---
  phi.current <- transform_fwd(theta.init)
  theta.current <- transform_inv(phi.current)
  logpost.current <- logPost(theta.current) + logJacobian(theta.current)

  chain <- matrix(NA, nrow=n.iter, ncol=p)
  acc <- 0

  # --- MCMC ---
  for(i in 1:n.iter){
    phi.prop <- mvrnorm(1, mu=phi.current, Sigma=Sigma)
    theta.prop <- transform_inv(phi.prop)
    logpost.prop <- logPost(theta.prop) + logJacobian(theta.prop)

    log.alpha <- logpost.prop - logpost.current
    if(log(runif(1)) < log.alpha){
      phi.current <- phi.prop
      theta.current <- theta.prop
      logpost.current <- logpost.prop
      acc <- acc + 1
    }
    chain[i,] <- theta.current
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
