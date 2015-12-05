#' @include psrf.R run_fixed_mcmc.R
NULL

#' @title Function \code{run_gelman_mcmc}
#' @description Accepts RNA-seq count data and group structure as arguments and 
#' returns a set of heterosis probabilities. 
#' @details Multiple MCMC chain is run until Gelman potential scale reduction factors
#' show no evidence of a lack of convergence.
#' @export
#' @return a \code{Chain} object
#' @param chain a \code{Chain} object
run_gelman_mcmc = function(chain){
  stopifnot(chain@nchains > 1)
  tol = 1.1

  if(chain@verbose){
    print(paste("Using Gelman-Rubin potential scale reduction factors on", chain@nchains, "separate chains to assess convergence."))
    print("Running pilot chain first to get dispersed starting values relative to the joint posterior.")
  }
 
  pilot = run_fixed_mcmc(chain)
  chain_list = list()
  chain_list[[1]] = pilot

  for(i in 2:chain@nchains){
    if(chain@verbose) print(paste0("Running dispersed chain ", i - 1, " of ", chain@nchains - 1, "."))
    chain_list[[i]] = run_fixed_mcmc(disperse_starts(pilot))
  }

  pilot@psrf = calc_gelman(chain_list)
  p = psrf(pilot, important = T, sort = T, threshold = 0)

  if(chain@verbose) {
    print("Summary of important Gelman factors:")
    print(summary(p))
    print("Highest 10:")   
    print(p[1:10])
    low = signif(mean(p < tol)* 100, 5)
    print(paste0(low, "% are less than ", tol, " (threshold)."))
    high = sum(p > tol)
    print(paste0(high, " are greater than ", tol, "."))
  }

  pilot
}
