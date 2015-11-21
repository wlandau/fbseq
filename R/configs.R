#' @include special_beta_priors.R parameters.R
NULL

#' @title Class \code{Configs}
#' @description A collection of MCMC control parameters.
#' @exportClass Configs
#'
#' @slot diag convergence diagnostic to use. Can be "gelman" or "none".
#' @slot ess Minimum effective sample size for all parameters
#' @slot max_attempts_diag Maximum number of retries for assessing convergence.
#' Can be set to Inf to run indefinitely.
#' @slot max_attempts_ess Maximum number of retries for generating enough effective samples.
#' Can be set to Inf to run indefinitely.
#' @slot nchains_diag number of independent chains to run (including this one) to use
#' convergence diagnostics that require multiple chains. 
#' @slot psrf_tol upper threshold for Gelman-Rubin potential scale reduction factors (if diag is "gelman")
#' 
#' @slot burnin MCMC burnin, the number of MCMC iterations to ignore at the beginning of each obj
#' @slot effects_update_beta values of l for which to update the beta_{l, g} parameters. For debugging only.
#' @slot genes_return Indices of genes whose parameter samples you want to return.
#' Applies to all gene-specific parameters except for the epsilons.
#' @slot genes_return_epsilon Indices of genes g for which epsilon_{n, g} is updated/returned.
#' @slot iterations Number of MCMC iterations (not including burnin or thinning)
#' @slot libraries_return Indices of RNA-seq libraries whose parameter samples you want to return.
#' @slot libraries_return_epsilon Indices of RNA-seq libraries n for which epsilon_{n, g} is updated/returned.
#' Applies to all library-specific parameters except for the epsilons.
#' @slot parameter_sets_return Character vector naming the variables whose MCMC samples 
#' you want to return
#' @slot parameter_sets_update Character vector naming the variables to calculate/update
#' during the MCMC.
#' @slot priors Names of the family of priors on the betas after integrating out the xi's. 
#' Can be any value returned by \code{special_beta_priors()}. All other values will default to the normal prior.
#' @slot thin MCMC thinning interval, number of iterations to skip in between iterations to return.
#' @slot verbose Number of times to print out progress during burnin and the actual MCMC.
#' If \code{verbose} > 0, then progress messages will also print during setup and cleanup.
setClass("Configs", 
  slots = list(
    diag = "character",
    ess = "numeric",
    max_attempts_diag = "numeric",
    max_attempts_ess = "numeric",
    nchains_diag = "numeric",
    psrf_tol = "numeric",

    burnin = "numeric",
    effects_update_beta = "numeric",
    genes_return = "numeric",
    genes_return_epsilon = "numeric",
    iterations = "numeric",
    libraries_return = "numeric",
    libraries_return_epsilon = "numeric",
    parameter_sets_return = "character",
    parameter_sets_update = "character",
    priors = "character",
    thin = "numeric",
    verbose = "numeric"
  ),
  prototype = list(
    diag = "gelman",
    ess = 1e2,
    max_attempts_diag = 1,
    max_attempts_ess = 0,
    nchains_diag = 4,
    psrf_tol = 1.1,

    burnin = 1e4,
    genes_return = numeric(0),
    genes_return_epsilon = numeric(0),
    iterations = 1e3,
    libraries_return = numeric(0),
    libraries_return_epsilon = numeric(0),
    parameter_sets_return = parameters(),
    parameter_sets_update = parameters(),
    priors = "normal",
    thin = 10,
    verbose = 5
  )
)

#' @title Constructor for class \code{Configs}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' Elements passed with \code{...} must be named. For example, \code{Configs(diag = "gelman")}.
#' @export
#' @return a \code{Configs} object
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... optional slot values.
Configs = function(obj = NULL, ...){
  configs = new("Configs", ...)

  if(class(obj) == "list") {
    for(n in intersect(slotNames(configs), names(obj)))
       slot(configs, n) = as(obj[[n]], class(slot(configs, n)))

    for(n in c("parameter_sets_return", "parameter_sets_update"))
      if(n %in% names(obj)){
        if(is.character(obj[[n]])){
          slot(configs, n) = as(obj[[n]], class(slot(configs, n)))
        } else {
          slot(configs, n) = as(names(obj[[n]])[as.logical(obj[[n]])], class(slot(configs, n)))
        }
      }
  } else if(class(obj) == "Chain"){
    subtract = c("parameter_sets_return", "parameter_sets_update", "priors")
    for(n in setdiff(intersect(slotNames(configs), slotNames(obj)), subtract))
      slot(configs, n) = as(slot(obj, n), class(slot(configs, n)))

    for(n in c("parameter_sets_return", "parameter_sets_update"))
      slot(configs, n) = as(names(slot(obj, n))[as.logical(slot(obj, n))], class(slot(configs, n)))

    configs@priors = as(ifelse(obj@priors > 0, special_beta_priors()[obj@priors], "normal"), class(configs@priors))
  }

  if(!any(configs@priors %in% special_beta_priors())){
    configs@parameter_sets_return = setdiff(configs@parameter_sets_return, "xi")
    configs@parameter_sets_update = setdiff(configs@parameter_sets_update, "xi")
  }
  
  configs
}