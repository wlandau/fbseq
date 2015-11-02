#' @include flatten.R
NULL

#' @title Function \code{credible_intervals}
#' @description Extracts posterior means and estimated credible intervals
#' from a \code{Chain} object. 
#'
#' @export
#' @return a data frame of estimates: posterior means and approximate credible intervals.
#' Parameters not updeted in the MCMC are excluded.
#' @param chain a \code{Chain} object
#' @param level level of the credible intervals from 0 to 1
credible_intervals = function(chain, level = 0.95){
  normals = c("beta", "epsilon", "theta")
  gammas = c("gamma", "nu", "omegaSquared", "rho", "sigmaSquared", "tau", "xi")

  Mean = flatten_post(chain)
  MeanSq = flatten_post(chain, square = T)
  Sd =  sqrt(chain@iterations*(MeanSq - Mean^2)/(chain@iterations - 1))

  d = data.frame(mean = Mean, sd = Sd, lower = NA, upper = NA)
  p = 1 - (1 - level)/2

  for(v in normals){
    n = grep(v, rownames(d))
    s = d[n,]
    s$lower = qnorm(1 - p, mean = s$mean, sd = s$sd)
    s$upper = qnorm(p, mean = s$mean, sd = s$sd)
    d[n,] = s
  }

  for(v in gammas){
    n = grep(v, rownames(d))
    s = d[n,]
    shape = s$mean^2/s$sd^2
    rate = s$mean/s$sd^2
    s$lower = qgamma(1 - p, shape = shape, rate = rate)
    s$upper = qgamma(p, shape = shape, rate = rate)
    d[n,] = s
  }

  d["nu", "upper"] = pmin(d["nu", "upper"], chain@d)
  d[grep("omegaSquared", rownames(d)), "upper"] = pmin(d[grep("omegaSquared", rownames(d)), "upper"], chain@w^2)
  d[grep("sigmaSquared", rownames(d)), "upper"] = pmin(d[grep("sigmaSquared", rownames(d)), "upper"], chain@s^2)

  out = data.frame(mean = d$mean, lower = d$lower, upper = d$upper)
  colnames(out) = c("mean", paste0(c("lower", "upper"), "_ci_", level))
  rownames(out) = names(Mean)
  out
}
