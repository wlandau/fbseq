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
#' @param backend Backend for parallel computing. Options are
#' "CUDA" and "OpenMP".
#' @param processes number of CPU processes to fork for 
#' the additional chains. This argument is automatically reset to 1 for non-serial
#' backends because additional \code{parallel::mclapply} processes interfere
#' with other modes of parallelism.
#' For some other backends, chains will be distributed accross processes.
#' @param threads Number of threads for the OpenMP implementation.
fbseq = function(chain, additional_chains = 3, backend = "CUDA", processes = 1, threads = 1){
  t = proc.time()
  processes = check_fbseq_input(backend, processes, threads)
  if(chain@verbose & additional_chains > 0) cat("Running pilot chain.\n")
  pilot = single_mcmc(chain, backend = backend, threads = threads)
  if(additional_chains < 1){
    out = pilot
  } else {
    if(chain@verbose){
      p = ifelse(processes > 1, 
        paste0("distributed over ", processes, " parallel processes.\n"), 
        "in sequence.\n")
      cat("Running", additional_chains, "additional dispersed chains", p)
    }
    out = c(pilot, mclapply(1:additional_chains + 1, function(i){
      dis = disperse_starts(pilot)
      if(chain@verbose) 
        cat("Running additional chain ", i - 1, " of ", additional_chains, ".\n", sep = "")
      single_mcmc(dis, backend = backend, threads = threads)
    }, mc.cores = processes))
  }
  attr(out, "runtime") = proc.time() - t
  return(out)
}
