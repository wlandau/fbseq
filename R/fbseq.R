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
#' @param backend defaults to "CUDA" (from package \code{fbseqCUDA}). 
#' Other options include "serial" (from package \code{fbseqSerial}), which
#' does not use any parallel computing.
#' @param processes number of CPU processes to fork for 
#' the additional chains. This argument is automatically reset to 1 for the CUDA 
#' backend because parallel processes can interfere with CUDA contexts.
#' For some other backends, chains will be distributed accross processes.
fbseq = function(chain, additional_chains = 3, backend = "CUDA", processes = 1){
  if(chain@verbose & additional_chains > 0) print(paste0("Running pilot chain with ", backend, " backend."))
  if(backend == "CUDA") processes = 1
  pilot = single_mcmc(chain, backend = backend)
  if(additional_chains < 1){
    return(pilot)
  } else {
    out = c(pilot, mclapply(1:additional_chains + 1, function(i){
      dis = disperse_starts(pilot)
      if(chain@verbose) print(paste0(paste0("Running additional chain ", i - 1, " of ", additional_chains, " with ", backend, " backend.")))
      single_mcmc(dis, backend = backend)
    }, mc.cores = processes))
    return(out)
  }
}
