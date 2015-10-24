#' @include cuda_usage.R
NULL

#' @title Function \code{run_fixed_mcmc}
#' @description Runs a single MCMC chain for a fixed number of iterations.
#' 
#' @export
#' @return a \code{Chain} object with updated parameter values. You can feed
#' this object back into another call to \code{run_fixed_mcmc()} to continue the chain
#' from where you left off.
#'
#' @param chain object of type \code{Chain}. This \code{chain} argument could be
#' a newly created \code{Chain} object from \code{Chain(...)}. Alternatively,
#' if \code{chain} is the output from a previous call to \code{run_mcmc(...)},
#' then the function will continue the MCMC from where it left off.
run_fixed_mcmc = function(chain){
  check_heterosisCUDA()
  chain@seeds = sample.int(1e3 * chain@N * chain@G, chain@N * chain@G)
  heterosisCUDA::heterosisCUDA(chain)
}
