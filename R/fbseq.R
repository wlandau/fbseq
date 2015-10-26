#' @include effectiveSampleSize.R run_gelman_mcmc.R
NULL

#' @title Function \code{fbseq}
#' @description Top-level function. Runs an MCMC with the runtime
#' parameters specified in the \code{Chain} object
#' 
#' @export
#' @return a \code{Chain} object with updated parameter values. You can feed
#' this object back into another call to \code{run_mcmc()} to continue the chain
#' from where you left off.
#'
#' @param chain object of type \code{Chain}. See the package vignette for details.
fbseq = function(chain){
  if(chain@diag == "gelman"){
    chain = run_gelman_mcmc(chain)
  } else {
    chain = run_fixed_mcmc(chain)
  }
  if(chain@ess >= 1) chain = effectiveSampleSize(chain)
  chain
}
