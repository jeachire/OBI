# OBI <img src="https://www.r-project.org/logo/Rlogo.png" width="40" align="right">

An R package for **Objective Bayesian Inference**.

## Installation

You can install the development version of **OBI** from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("jeachire/OBI")
```

Here is a minimal example using the main function `fitBayes()`:

```r
library(OBI)

# Simulate censored Weibull data
set.seed(123)
x <- rweibull(100, shape = 0.8, scale = 2)     # lifetimes
censor <- runif(100, 0, 3)                     # censoring times
delta <- as.numeric(x <= censor)               # censoring indicator
x_cens <- pmin(x, censor)                      # observed data

# Fit Bayesian Weibull model with Jeffreys prior
fit <- fitBayes(
  x = x_cens,
  dist = "weibull",
  delta = delta,
  prior = "Jeffreys",
  method = "MCMC",
  n.iter = 5000,
  burnin = 1000,
  thin = 2
)

# Posterior summary
fit$summary

# Trace plot for parameters
plot(fit$posterior[,1], type = "l", main = "Shape parameter")
plot(fit$posterior[,2], type = "l", main = "Scale parameter")
````
