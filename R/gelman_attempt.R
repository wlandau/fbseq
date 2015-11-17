#' @include run_fixed_mcmc.R
NULL

#' @title Function \code{gelman_attempt}
#' @description Run 4 chains together and compute and report Gelman factors. 
#' @export
#' @return a list of \code{Chain} objects
#' @param chain_list a list of \code{Chain} objects
#' @param pattern character string: grep pattern for the parameters of interest
gelman_attempt = function(chain_list, pattern){
  if(chain_list[[1]]@verbose) print(paste0("Generating dispersed starting values for the other ", chain_list[[1]]@nchains_diag - 1," chains using pilot chain."))
  for(i in 2:chain_list[[1]]@nchains_diag) chain_list[[i]] = disperse_starts(chain_list[[1]])

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