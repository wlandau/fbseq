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

#' @title Function \code{gelman_attempt}
#' @description Run 4 chains together and compute and report Gelman factors. 
#' @export
#' @return a list of \code{Chain} objects
#' @param chain_list a list of \code{Chain} objects
#' @param pattern character string: grep pattern for the parameters of interest
gelman_attempt = function(chain_list, pattern){
  chain_list = lapply(chain_list, run_fixed_mcmc)
  all_psrf = calc_gelman(chain_list)
  all_psrf[!is.finite(all_psrf)] = 0
  for(i in 1:chain_list[[1]]@nchains_diag) chain_list[[i]]@psrf = all_psrf
  psrf = all_psrf[grepl(pattern, names(all_psrf))]
  tol = chain_list[[1]]@psrf_tol

  if(chain_list[[1]]@verbose) {
    print("Summary of Gelman factors for betas involved in contrasts and hyperparameters:")
    print(summary(psrf))
    print("Highest 10:")   
    print(sort(psrf, decreasing = T)[1:10])
    low = signif(mean(psrf < tol)* 100, 5)
    print(paste0(low, "% are less than ", tol, " (threshold)."))
    high = sum(psrf > tol)
    print(paste0(high, " are greater than ", tol, "."))
  }

  chain_list
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

  if(chain@verbose) print(paste0("Generating dispersed starting values for the other ", chain@nchains_diag - 1," chains using pilot chain."))
  chain_list = list(pilot_chain)
  for(i in 2:chain@nchains_diag) chain_list[[i]] = disperse_starts(pilot_chain)
 
  betas = which(apply(do.call(rbind, Scenario(chain)@contrasts), 2, function(x) any(x != 0)))
  pattern = paste(c(paste0("beta_", betas), "nu",  "omegaSquared", "sigmaSquared", "tau", "theta"), collapse = "|")

  attempt = 0
  while(attempt < chain@max_attempts_redisperse){
    attempt = attempt + 1
    if(chain@verbose) print(paste0("Attempt ", attempt, " of ", chain@max_attempts_redisperse,
      ": running ", chain@nchains_diag, " parallel chains. Regenerating dispersed starting values for the other chains relative to the pilot chain before every attempt."))
    chain_list = gelman_attempt(chain_list, pattern)
    psrf = chain_list[[1]]@psrf
    if(all(psrf[grepl(pattern, names(psrf))] < chain@psrf_tol, na.rm = T)) break
    if(chain@verbose) print(paste0("Regenerating dispersed starting values for the other ", chain@nchains_diag - 1," chains relative to the current pilot chain."))
    for(i in 2:chain@nchains_diag) chain_list[[i]] = disperse_starts(chain_list[[1]])
  }
  chain_list[[1]]@attempts_redisperse = as.integer(attempt)
  if(attempt == chain@max_attempts_redisperse) warning(paste("In run_gelman_mcmc(), chain@max_attempts_redisperse =", chain@max_attempts_redisperse, "reached."))

  chain_list = lapply(chain_list, function(ch){ch@burnin = as.integer(0); ch})
  attempt = 0

  while(attempt < chain@max_attempts_diag){
    attempt = attempt + 1
    if(chain@verbose) print(paste0("Attempt ", attempt, " of ", chain@max_attempts_diag,
      ": running ", chain@nchains_diag, " parallel chains without redispersing starting values."))
    chain_list = gelman_attempt(chain_list, pattern)
    psrf = chain_list[[1]]@psrf
    if(all(psrf[grepl(pattern, names(psrf))] < chain@psrf_tol, na.rm = T)) break
  }
  chain_list[[1]]@attempts_diag = as.integer(attempt)
  if(attempt == chain@max_attempts_diag) warning(paste("In run_gelman_mcmc(), chain@max_attempts_diag =", chain@max_attempts_diag, "reached."))

  chain = chain_list[[1]]
  chain@burnin = as.integer(burnin)
  chain
}
