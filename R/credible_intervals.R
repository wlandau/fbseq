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
  Mean = flatten_post(chain)
  MeanSq = flatten_post(chain, square = T)
  Sd =  sqrt(chain@iterations*(MeanSq - Mean^2)/(chain@iterations - 1))
  p = 1 - (1 - level)/2

  lowerCInorm = qnorm(1 - p, mean = Mean, sd = Sd)
  upperCInorm = qnorm(p, mean = Mean, sd = Sd)

  shapeGamma = abs(Mean)^2/Sd^2
  rateGamma = abs(Mean)/Sd^2
  lowerCIgamma = qgamma(1 - p, shape = shapeGamma, rate = rateGamma)
  upperCIgamma = qgamma(p, shape = shapeGamma, rate = rateGamma)

  lowerCI = lowerCInorm
  upperCI = upperCInorm
  
  useGamma = NULL
  for(parm in c("gamma", "nu", "omegaSquared", "rho", "sigmaSquared", "tau", "xi"))
    useGamma = c(useGamma, grep(parm, names(Mean)))

  lowerCI[useGamma] = lowerCIgamma[useGamma]
  upperCI[useGamma] = upperCIgamma[useGamma]

  if(is.finite(upperCI["nu"])) upperCI["nu"] = min(upperCI["nu"], chain@d)
  if(is.finite(upperCI["omegaSquared"])) upperCI["omegaSquared"] = min(upperCI["omegaSquared"], chain@w^2)
  if(is.finite(upperCI["sigmaSquared"])) upperCI["sigmaSquared"] = min(upperCI["sigmaSquared"], chain@s^2)

  lowerCI[grep("delta", names(Mean))] = 0
  upperCI[grep("delta", names(Mean))] = 1

  lowerCI[grep("pi", names(Mean))] = max(0, lowerCI[grep("pi", names(Mean))])
  upperCI[grep("pi", names(Mean))] = min(1, upperCI[grep("pi", names(Mean))])

  out = data.frame(Mean, lowerCI, upperCI)
  colnames(out) = c("mean", paste0(c("lower", "upper"), "_ci_", level))
  out
}
