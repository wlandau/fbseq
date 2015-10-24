#' @include alternate_priors.R parameters.R
NULL

#' @title Class \code{Configs}
#' @description A collection of MCMC control parameters.
#' @exportClass Configs
#' 
#' @slot diag convergence diagnostic to use. Can be "gelman" or "none".
#' @slot ess Minimum effective sample size for all parameters
#' @slot psrf_tol upper threshold for Gelman-Rubin potential scale reduction factors (if diag is "gelman")
#' @slot max_attempts Maximum number of retries for assessing convergence and generating enough effective samples.
#' Can be set to Inf to run indefinitely.
#' @slot nchains_diag number of independent chains to run (including this one) to use
#' convergence diagnostics that require multiple chains. 
#' 
#' @slot iterations Number of MCMC iterations (not including burnin or thinning)
#' @slot burnin MCMC burnin, the number of MCMC iterations to ignore at the beginning of each obj
#' @slot thin MCMC thinning interval, number of iterations to skip in between iterations to return.
#' @slot verbose Number of times to print out progress during burnin and the actual MCMC.
#' If \code{verbose} > 0, then progress messages will also print during setup and cleanup.
#' @slot mphtol Tolerance theshold for mid-parent heterosis (on a log base e scale)
#' 
#' @slot returns Character vector naming the variables whose MCMC samples 
#' you want to return
#' @slot updates Character vector naming the variables to calculate/update
#' during the MCMC.
#' 
#' @slot samples_return Indices of RNA-seq samples/libraries whose parameter samples you want to return.
#' Applies to all library-specific parameters except for the epsilons.
#' @slot features_return Indices of features/genes whose parameter samples you want to return.
#' Applies to all gene-specific parameters except for the epsilons.
#'
#' @slot samples_return_eps Indices of RNA-seq samples/libraries n for which epsilon_{n, g} is updated/returned.
#' @slot features_return_eps Indices of features/genes g for which epsilon_{n, g} is updated/returned.
#'
#' @slot phiPrior Name of the family of priors on the phi_g's after integrating out the xi_phi's. 
#' Can be any value returned by \code{alternate_priors()}. All other values will default to the normal prior.
#' @slot alpPrior Name of the family of priors on the alp_g's after integrating out the xi_alp's. 
#' Can be any value returned by \code{alternate_priors()}. All other values will default to the normal prior.
#' @slot delPrior Name of the family of priors on the del_g's after integrating out the xi_del's. 
#' Can be any value returned by \code{alternate_priors()}. All other values will default to the normal prior.
setClass("Configs", 
  slots = list(
    diag = "character",
    ess = "numeric",
    psrf_tol = "numeric",
    max_attempts = "numeric",
    nchains_diag = "numeric",
    iterations = "numeric",
    burnin = "numeric",
    thin = "numeric",
    verbose = "numeric",
    mphtol = "numeric",
    returns = "character",
    updates = "character",
    samples_return = "numeric",
    features_return = "numeric",
    samples_return_eps = "numeric",
    features_return_eps = "numeric",
    phiPrior = "character",
    alpPrior = "character",
    delPrior = "character"
  ),
  prototype = list(
    diag = "gelman",
    ess = 1e2,
    psrf_tol = 1.1,
    max_attempts = 10,
    nchains_diag = 4,
    iterations = 1e3,
    burnin = 1e4,
    thin = 1e1,
    verbose = 1e1,
    mphtol = 5e-1,
    returns = setdiff(parameters(), "tauGam"),
    updates = setdiff(parameters(), "tauGam"),
    phiPrior = "normal",
    alpPrior = "normal",
    delPrior = "normal"
  )
)

#' @title Constructor for class \code{Configs}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' @export
#' @return a \code{Configs} object
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... optional slot values.
Configs = function(obj = NULL, ...){
  configs = new("Configs", ...)

  if(class(obj) == "list") {
    for(n in intersect(slotNames(configs), names(obj)))
       slot(configs, n) = obj[[n]]

    if("M" %in% names(obj))
      configs@iterations = obj$M

    for(n in c("returns", "updates"))
      if(n %in% names(obj))
        slot(configs, n) = names(obj[[n]])[as.logical(obj[[n]])]
  } else if(class(obj) == "Chain"){
    priors = paste0(c("phi", "alp", "del"), "Prior")
    subtract = c("returns", "updates", priors)
    for(n in setdiff(intersect(slotNames(configs), slotNames(obj)), subtract))
      slot(configs, n) = slot(obj, n)

    for(n in c("returns", "updates"))
      slot(configs, n) = names(slot(obj, n))[as.logical(slot(obj, n))]

    for(n in priors)
      if(slot(obj, n)){
        slot(configs, n) = alternate_priors()[slot(obj, n)]
      } else {
        slot(configs, n) = "normal"
      }

    configs@iterations = obj@M
  }

  for(p in c("phi", "alp", "del")){
    cap = paste0(toupper(substring(p, 1, 1)), substring(p, 2),collapse="")
    xi = paste0("xi", cap)
    if(!(slot(configs, paste0(p, "Prior")) %in% alternate_priors())){
      configs@returns = setdiff(configs@returns, xi)
      configs@updates = setdiff(configs@updates, xi)
    }
  }

  return(configs)
}