#' @include effectiveSampleSize.R flatten.R run_fixed_mcmc.R
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

  if(chain@burnin < 1e4) warning("small burnin. Chains generated with disperse_starts() may be numerically unstable.")
  pilot_chain = run_fixed_mcmc(chain)

  if(chain@verbose) print(paste0("Generating dispersed starting values for the other ", chain@nchains_diag - 1," chains using pilot chain."))
  chain_list = list(pilot_chain)
  for(i in 2:chain@nchains_diag) chain_list[[i]] = disperse_starts(pilot_chain)
 
  attempt = 0
  pattern = paste(c("beta", "nu", "omega", "tau", "theta"), collapse = "|")

  while(attempt < chain@max_attempts){
    attempt = attempt + 1
    if(chain@verbose) print(paste0("Attempt ", attempt, " of ", chain@max_attempts,
      ": running ", chain@nchains_diag, " parallel chains."))
    chain_list = lapply(chain_list, run_fixed_mcmc)

    all_psrf = calc_gelman(chain_list)
    all_psrf[!is.finite(all_psrf)] = 0
    for(i in 1:chain@nchains_diag) chain_list[[i]]@psrf = all_psrf
    psrf = all_psrf[grepl(pattern, names(all_psrf))]

    if(chain@verbose) {
      print("Summary of Gelman factors for hyperparameters and betas:")
      print(summary(psrf))
      print("Highest 10 Gelman factors:")   
      print(sort(psrf, decreasing = T)[1:10])
      low = signif(mean(psrf < chain@psrf_tol)* 100, 5)
      print(paste0(low, "% of Gelman factors are less than ", chain@psrf_tol, " (threshold)."))
      high = sum(psrf > chain@psrf_tol)
      print(paste0(high, " Gelman factors are greater than ", chain@psrf_tol, "."))
    }

    if(any(psrf > chain@psrf_tol, na.rm = T)){
      chain_list = lapply(chain_list, function(chain){chain@burnin = as.integer(0); chain})
    } else {
      break
    }
  }

  chain = chain_list[[1]]
  iterations = chain@iterations
  for(i in 2:chain@nchains_diag) chain = concatenate(chain_list[[i]], chain)
  if(attempt == chain@max_attempts) warning(paste("In run_gelman_mcmc(), chain@max_attempts =", chain@max_attempts, "reached."))
  chain
}
