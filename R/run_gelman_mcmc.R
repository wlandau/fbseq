#' @include run_fixed_mcmc.R
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
  stopifnot(chain@nchains_diag > 1)

  betas = which(apply(do.call(rbind, Scenario(chain)@contrasts), 2, function(x) any(x != 0)))
  pattern = paste(c(paste0("beta_", betas), "nu",  "omegaSquared", "sigmaSquared", "tau", "theta"), collapse = "|")

  if(chain@verbose){
    print(paste("Using Gelman-Rubin potential scale reduction factors on", chain@nchains_diag, "parallel chains to assess convergence."))
    print("Running pilot chain first to get dispersed starting values relative to the joint posterior.")
  }
 
  chain_list = list(run_fixed_mcmc(chain))

  attempt = 0
  while(attempt < chain@max_attempts_diag){
    attempt = attempt + 1
    if(chain_list[[1]]@verbose) print(paste0("Attempt ", attempt, " of ", chain_list[[1]]@max_attempts_diag,
    ": running ", chain_list[[1]]@nchains_diag, " parallel chains."))

    chain_list = gelman_attempt(chain_list, pattern)
    if(all(chain_list[[1]]@psrf[grepl(pattern, names(chain_list[[1]]@psrf))] < chain@psrf_tol, na.rm = T)) break
    if(attempt < chain@max_attempts_diag)
      chain_list = lapply(chain_list, function(ch){
        ch@burnin = as.integer(2*ch@burnin)
        ch@thin = as.integer(2*ch@thin)
        ch
      })
  }

  chain_list[[1]]@attempts_diag = as.integer(attempt)
  if(attempt == chain@max_attempts_diag) warning(paste("In run_gelman_mcmc(), chain@max_attempts_diag =", chain@max_attempts_diag, "reached."))
  chain_list[[1]]
}
