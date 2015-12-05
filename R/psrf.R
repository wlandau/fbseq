#' @title Function \code{psrf}
#' @description Get the (important) Gelman-Rubin potential scale reduction factors (psrf) from a \code{Chain}
#' object.
#' @export
#' @return a vector of important potential scale reduction factors
#' @param chain a \code{Chain} object
#' @param important If \code{TRUE}, only the important psrf will be returned. The 
#' important psrf correspond to the hyperparameters (thetas, sigmas, nu, and tau)
#' plus the betas used in estimating the posterior probabilities specified in \code{Scenario(chain)}.
#' @param sort sort the psrf in decreasing order
#' @param threshold If only the psrf above \code{threshold} will be returned. To return all psrf,
#' set \code{threshold} to 0.
psrf = function(chain, important = T, sort = T, threshold = 1.1){
  psrf = chain@psrf
  betas = which(apply(do.call(rbind, Scenario(chain)@contrasts), 2, function(x) any(x != 0)))
  pattern = paste(c(paste0("beta_", betas), "nu", "sigmaSquared", "tau", "theta"), collapse = "|")
  if(important) psrf = psrf[grepl(pattern, names(psrf))]
  if(sort) psrf = sort(psrf, decreasing = T)
  psrf[psrf > threshold]
}
