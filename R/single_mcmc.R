#' @include check_backend.R cuda_usage.R
NULL

#' @title Function \code{single_mcmc}
#' @description Runs a single MCMC chain for a fixed number of iterations.
#' 
#' @export
#' @return a \code{Chain} object with updated parameter values. You can feed
#' this object back into another call to \code{single_mcmc()} to continue the chain
#' from where you left off.
#'
#' @param chain object of type \code{Chain}. This \code{chain} argument could be
#' a newly created \code{Chain} object from \code{Chain(...)}. Alternatively,
#' if \code{chain} is the output from a previous call to \code{run_mcmc(...)},
#' then the function will continue the MCMC from where it left off.
#' @param backend defaults to "CUDA". The other option is "OpenMP",
#' which uses \code{OpenMP} threads to parallelize within Markov chains
#' if \code{OpenMP} is available.
#' @param threads number of threads for the OpenMP backend.
single_mcmc = function(chain, backend = "CUDA", threads = 1){
  chain@seeds = as.integer(sample.int(1e3 * chain@N * chain@G, chain@N * chain@G))
  if(chain@verbose) cat("Using", backend, "backend.\n")
  if(backend == "CUDA"){
    check_backend("fbseqCUDA")
    fbseqCUDA::fbseqCUDA(chain)
  } else if(backend == "OpenMP"){
    check_backend("fbseqOpenMP")
    fbseqOpenMP::fbseqOpenMP(chain, threads = threads)
  } else {
    stop("illegal backend.")
  }
}
