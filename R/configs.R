#' @include parameters.R special_beta_priors.R
NULL

#' @title Class \code{Configs}
#' @description Set individual MCMC control parameters 
#' (slots listed in \code{help("Configs-class")}).
#' @seealso \code{help("Configs-class")}
#' @exportClass Configs
#'
#' @slot burnin MCMC burnin, the number of MCMC iterations to ignore at the beginning of each obj
#' @slot effects_update_beta values of l for which to update the beta_{l, g} parameters. For debugging only.
#' @slot genes_return Indices of genes whose parameter samples you want to return.
#' Applies to all gene-specific parameters except for the epsilons.
#' @slot genes_return_epsilon Indices of genes g for which epsilon_{n, g} is updated/returned.
#' @slot iterations Number of MCMC iterations after burnin for which selected parameter samples are kept.
#' Total MCMC iterations = burnin + thin * "iterations", and the whole "thin * iterations" portion
#' is used to calculate posterior means, mean squares, and probabilities.
#' @slot libraries_return Indices of RNA-seq libraries whose parameter samples you want to return.
#' Currently moot because there are no library-specific parameters other than the epsilons, but that
#' could change in future versions of the package.
#' @slot libraries_return_epsilon Indices of RNA-seq libraries n for which epsilon_{n, g} is updated/returned.
#' Applies to all library-specific parameters except for the epsilons.
#' @slot parameter_sets_return Character vector naming the variables whose MCMC samples 
#' you want to return
#' @slot parameter_sets_update Character vector naming the variables to calculate/update
#' during the MCMC.
#' @slot priors Character vector. Names of the family of priors on the betas after integrating out the xi's. 
#' Can be any value returned by \code{special_beta_priors()}. All other values will default to the normal prior.
#' @slot samplers character string indicating the sampling algorithm
#' @slot thin MCMC thinning interval. \code{thin = 1} means parameter samples will be saved for every iterations
#' after burnin. \code{thin = 10} means parameter samples will be saved every 10th iteration after burnin.
#' Total MCMC iterations = burnin + thin * "iterations", and the whole "thin * iterations" portion
#' is used to calculate posterior means, mean squares, and probabilities.
#' @slot verbose Number of times to print out progress during burnin and the actual MCMC.
#' If \code{verbose} > 0, then progress messages will also print during setup and cleanup.
setClass("Configs", 
  slots = list(
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
    samplers = "character",
    thin = "numeric",
    verbose = "numeric"
  ),
  prototype = list(
    burnin = 1e5,
    genes_return = numeric(0),
    genes_return_epsilon = numeric(0),
    iterations = 5e3,
    libraries_return = numeric(0),
    libraries_return_epsilon = numeric(0),
    parameter_sets_return = parameters(),
    parameter_sets_update = parameters(),
    priors = "normal",
    samplers = "default",
    thin = 20,
    verbose = 5
  )
)

#' @title Constructor for class \code{Configs}
#' @description Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
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

  configs@thin = max(1, configs@thin)

  if(!any(configs@priors %in% special_beta_priors())){
    configs@parameter_sets_return = setdiff(configs@parameter_sets_return, "xi")
    configs@parameter_sets_update = setdiff(configs@parameter_sets_update, "xi")
  }
  
  configs
}