#' @include concatenate.R run_fixed_mcmc.R
NULL

#' @title Function \code{effectiveSampleSize}
#' @description Reruns an already-run \code{Chain} object until 
#' at least \code{chain@@M} effective samples are obtained for each returned
#' parameter. 
#' @export
#' @return a \code{Chain} object
#' @param chain a \code{Chain} object
effectiveSampleSize = function(chain){
  i = 0
  iterations = chain@M
  while(i < chain@max_attempts){
    i = i + 1
    flat = flatten(chain)

    if(any(dim(flat) < 1)){
      if(chain@verbose) print("No returned MCMC samples. Returning from effectiveSampleSize().")
      return(chain)
    }    

    ess = effectiveSize(mcmc(flat))
    hyper = intersect(hyperparameters(), names(ess))

    if(chain@verbose) {
      print(paste0("Attempt ", i, " of ", chain@max_attempts, ": trying to obtain ", chain@ess, " effective samples for every parameter with returned samples."))
      print("Summary of effective sample sizes of returned parameters:")
      print(summary(ess))
      print("Effective sample sizes of returned hyperparameters:")      
      print(ess[hyper])
    }

    min_ess = round(min(ess), 3)
    which_min_ess = names(ess)[which.min(ess)]
    if(min_ess >= chain@ess) return(chain)

    if(chain@verbose) print(paste0("Only ", min_ess, " effective samples obtained for ", which_min_ess, ". Continuing until ", chain@ess, " effective samples reached."))

    configs = Configs(chain)
    configs@burnin = 0
    configs@iterations = iterations
    new_chain = Chain(matrix(chain@counts, ncol = chain@N), chain@group, configs, Starts(chain))
    new_chain = run_fixed_mcmc(new_chain)
    chain = concatenate(chain, new_chain)
  }

  warning(paste("In effectiveSampleSize(), chain@max_attempts =", chain@max_attempts, "reached."))
  chain
}
