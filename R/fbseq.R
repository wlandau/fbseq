#' @include single_mcmc.R
NULL

#' @title Function \code{fbseq}
#' @description Top-level function. Runs an MCMC with the runtime
#' parameters specified in the \code{Chain} object
#' @export
#' @return a \code{Chain} object or a list of \code{Chain} objects, depending on 
#' the value of \code{additional_chains}.
#' @param chain object of type \code{Chain}. See the package vignette for details.
#' @param additional_chains If > 0, \code{additional_chains} additional \code{Chain} objects will be created from the 
#' original \code{chain} argument and run. This is useful for calculating convergence diagnostics
#' such as Gelman-Rubin potential scale reduction factors. The additional chains will be run after
#' the first chain and have starting values overdispersed to the full joint posterior distribution, as estimated
#' by the results of the first chain.
fbseq = function(chain, additional_chains = 3){
  if(chain@verbose & additional_chains > 0) print("Running pilot chain.")
  pilot = single_mcmc(chain)
  if(additional_chains < 1){
    return(pilot)
  } else {
    out = list(pilot)
    for(i in 1:additional_chains + 1){
      print(paste0("Running additional chain ", i, " of ", additional_chains, "."))
      out[[i]] = single_mcmc(disperse_starts(pilot))
    }
    return(out)
  }
}
