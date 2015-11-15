#' @include flatten.R run_fixed_mcmc.R
NULL

#' @title Function \code{calc_gelman}
#' @description calculate Gelman potential scale reduction factors on a list of parallel independent chains.
#' @export
#' @return Gelman potential scale reduction factors on a list of parallel independent chains.
#' @param chain_list list of \code{Chain} objects
calc_gelman = function(chain_list){
  stopifnot(length(chain_list) > 1)
  stopifnot(length(unique(sapply(chain_list, function(chain) chain@iterations))) == 1)

  Mean = sapply(chain_list, flatten_post)
  MeanSquare = sapply(chain_list, flatten_post, square = T)

  m = length(chain_list)
  n = chain_list[[1]]@iterations

  SjSq =  n*(MeanSquare - Mean^2)/(n - 1)
  W = rowMeans(SjSq)
  B = n*apply(Mean, 1, var)
  sqrt((n-1 + B/W)/n)
}

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
  if(chain@verbose){
    print(paste("Using Gelman-Rubin potential scale reduction factors on", chain@nchains_diag, "parallel chains to assess convergence."))
    print("Running pilot chain first to get dispersed starting values relative to the joint posterior.")
  }

  if(chain@burnin < 1e5) warning("small burnin. Chains generated with disperse_starts() may be slow to converge.")
  pilot_chain = run_fixed_mcmc(chain)
  burnin = pilot_chain@burnin
  pilot_chain@burnin = pilot_chain@burnin_diag

  if(chain@verbose) print(paste0("Generating dispersed starting values for the other ", chain@nchains_diag - 1," chains using pilot chain."))
  chain_list = list(pilot_chain)
  for(i in 2:chain@nchains_diag) chain_list[[i]] = disperse_starts(pilot_chain)
 
  attempt = 0
  betas = which(apply(do.call(rbind, Scenario(chain)@contrasts), 2, function(x) any(x != 0)))
  pattern = paste(c(paste0("beta_", betas), "nu",  "omegaSquared", "sigmaSquared", "tau", "theta"), collapse = "|")

  while(attempt < chain@max_attempts_diag){
    attempt = attempt + 1
    if(chain@verbose) print(paste0("Attempt ", attempt, " of ", chain@max_attempts_diag,
      ": running ", chain@nchains_diag, " parallel chains."))
    chain_list = lapply(chain_list, run_fixed_mcmc)

    all_psrf = calc_gelman(chain_list)
    all_psrf[!is.finite(all_psrf)] = 0
    for(i in 1:chain@nchains_diag) chain_list[[i]]@psrf = all_psrf
    psrf = all_psrf[grepl(pattern, names(all_psrf))]

    if(chain@verbose) {
      print("Summary of Gelman factors for betas involved in contrasts and hyperparameters:")
      print(summary(psrf))
      print("Highest 10:")   
      print(sort(psrf, decreasing = T)[1:10])
      low = signif(mean(psrf < chain@psrf_tol)* 100, 5)
      print(paste0(low, "% are less than ", chain@psrf_tol, " (threshold)."))
      high = sum(psrf > chain@psrf_tol)
      print(paste0(high, " are greater than ", chain@psrf_tol, "."))
    }

    if(all(psrf < chain@psrf_tol, na.rm = T)) break
  }

  chain = chain_list[[1]]
  if(attempt == chain@max_attempts_diag) warning(paste("In run_gelman_mcmc(), chain@max_attempts_diag =", chain@max_attempts_diag, "reached."))
  pilot_chain@burnin = as.integer(burnin)
  chain@attempts_diag = as.integer(attempt)
  chain
}
