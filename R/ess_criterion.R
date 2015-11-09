#' @include concatenate.R run_fixed_mcmc.R
NULL

#' @title Function \code{ess_criterion}
#' @description Reruns an already-run \code{Chain} object until 
#' at least \code{chain@@iterations} effective samples are obtained for each returned hyperparameter.
#' @export
#' @return a \code{Chain} object
#' @param chain a \code{Chain} object
ess_criterion = function(chain){
  i = 0
  iterations = chain@iterations
  while(i < chain@max_attempts){
    i = i + 1
    flat = flatten(chain)

    if(any(dim(flat) < 1)){
      if(chain@verbose) print("No returned MCMC samples. Returning from effectiveSampleSize().")
      return(chain)
    }    

    ess = effectiveSize(mcmc(flat))
    pattern = paste(c("nu", "omegaSquared", "sigmaSquared", "tau", "theta"), collapse = "|")
    ess = ess[grep(pattern, names(ess))]

    if(chain@verbose) {
      print(paste0("Attempt ", i, " of ", chain@max_attempts, ": trying to obtain ", chain@ess, " effective samples for every parameter with returned samples."))
      print("Summary of effective sample sizes of returned parameters:")
      print(summary(ess))
    }

    min_ess = round(min(ess), 3)
    which_min_ess = names(ess)[which.min(ess)]
    if(min_ess >= chain@ess) return(chain)

    if(chain@verbose) print(paste0("Only ", min_ess, " effective samples obtained for ", which_min_ess, ". Continuing until ", chain@ess, " effective samples reached."))

    configs = Configs(chain)
    configs@burnin = 0
    configs@iterations = iterations
    new_chain = Chain(Scenario(chain), configs, Starts(chain))
    new_chain = run_fixed_mcmc(new_chain)
    chain = concatenate(chain, new_chain)
  }

  warning(paste("In effectiveSampleSize(), chain@max_attempts =", chain@max_attempts, "reached."))
  chain
}
