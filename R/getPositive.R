#' @title Identify Positive Parameter Constraints
#' @description Returns a logical vector indicating which parameters of a given
#' distribution are restricted to the positive real line \eqn{(0, \infty)}.
#'
#' @param dist A character string specifying the distribution name.
#' Supported options include: \code{"normal"}, \code{"poisson"},
#' \code{"exponential"}, \code{"gamma"}, and \code{"weibull"}.
#'
#' @return A logical vector of length equal to the number of parameters
#' of the specified distribution. Each element indicates whether the
#' corresponding parameter is restricted to \eqn{(0, \infty)}
#' (\code{TRUE}) or not (\code{FALSE}).
#'
#' @details
#' - \code{"normal"}: parameters = mean, sd → \code{c(FALSE, TRUE)}
#' - \code{"poisson"}: parameter = lambda → \code{c(TRUE)}
#' - \code{"exponential"}: parameter = rate → \code{c(TRUE)}
#' - \code{"gamma"}: parameters = shape, rate → \code{c(TRUE, TRUE)}
#' - \code{"weibull"}: parameters = shape, scale → \code{c(TRUE, TRUE)}
#'
#' If an unsupported distribution is provided, the function throws an error.
#'
#' @examples
#' getPositive("normal")
#' getPositive("gamma")
#' getPositive("weibull")
#'
#' @export
getPositive <- function(dist) {
  if (dist == "normal") {
    # parameters: mean, sd
    return(c(FALSE, TRUE))
  } else if (dist == "poisson") {
    # parameter: lambda
    return(c(TRUE))
  } else if (dist == "exponential") {
    # parameter: rate
    return(c(TRUE))
  } else if (dist == "gamma") {
    # parameters: shape, rate
    return(c(TRUE, TRUE))
  } else if (dist == "weibull") {
    # parameters: shape, scale
    return(c(TRUE, TRUE))
  } else {
    stop("Unsupported distribution in getPositive().")
  }
}
