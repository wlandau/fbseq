#' @include alternate_priors.R generate_starts.R
NULL

#' @title Function \code{Chain}
#' @description Construct a \code{Chain} object, the 
#' one argument to \code{single_chain(...)}.
#' @details If \code{list} is set, all other arguments will be dropped.
#' @export
#' @return a \code{Chain} object ready to be passed to \code{single_chain()}
#'
#' @param data A data frame/matrix of RNA-seq read counts or list of slots. If a list,
#' then the function will return a \code{Chain} object with those slots. Note: \code{data}
#' should only be a list within internal functions of the package. It is not recommended that
#' the user assign \code{data} to be a list.
#' @param group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for parent 2,
#' and 3 for the hybrid.
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting values.
Chain = function(
  data, 
  group, 
  configs = Configs(),
  starts = Starts()
){
  chain = new("Chain")

  if(class(data) == "list"){
    for(n in intersect(names(data), slotNames(chain)))
       slot(chain, n) = data[[n]]
    return(chain)
  }

  counts = as.matrix(data)
  stopifnot(all(is.finite(counts)))
  chain = plug_in_chain(chain, configs = Configs(), starts = Starts())
  chain = plug_in_chain(chain, configs = configs, starts = starts)
  chain = fill_easy_gaps(chain, counts, group)
  chain = simple_starts(chain, counts, group)
  if(chain@returns["tauGam"] || chain@updates["tauGam"]) warning("TauGam should be set constant at 1. Sampling tauGam will prevent the MCMC from converging, and returning even a constant tauGam will muddle effective sample size diagnostics. Control tauGam with the returns, updates, returns_skip, and updates_skip slots in your Configs object. See the package vignettes for details.")
  chain
}

#' @title Function \code{plug_in_chain}
#' @description Plug the slots of \code{Configs} and \code{Starts}
#' objects into the analogous slots of a \code{Chain} object. Used
#' internally
#' @export
#' @return a \code{Chain} object 
#'
#' @param chain a \code{Chain} object
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting values.
plug_in_chain = function(chain, configs, starts){
  chain@M = as.integer(configs@iterations)
  priors = paste0(c("phi", "alp", "del"), "Prior")
  subtract = c("returns", "updates", priors)
  for(n in setdiff(intersect(slotNames(chain), slotNames(configs)), subtract))
    slot(chain, n) = as(slot(configs, n), class(slot(chain, n)))

  for(n in c("returns", "updates")){
    slot(chain, n) = as.integer(parameters() %in% slot(configs, n))
    names(slot(chain, n)) = parameters()
  }

  for(n in slotNames(chain)){
    if(grepl("Start", n))
      slot(chain, n) = as(slot(starts, gsub("Start", "", n)), class(slot(chain, n)))
    else if(n %in% slotNames(starts) && !(paste(n, "Start", sep = "") %in% slotNames(chain)))
      slot(chain, n) = as(slot(starts, n), class(slot(chain, n)))
  }
  
  for(n in priors)
    if(slot(configs, n) %in% alternate_priors()){
      slot(chain, n) = which(alternate_priors() == slot(configs, n))
    } else {
      slot(chain, n) = as.integer(0)
    }

  return(chain)
}

#' @title Function \code{fill_easy_gaps}
#' @description Set required slots in a \code{Chain} object.
#' Used internally.
#' @export
#' @return a \code{Chain} object 
#'
#' @param chain a \code{Chain} object
#' @param counts Matrix of RNA-seq read counts.
#' @param group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for parent 2,
#' and 3 for the hybrid.
fill_easy_gaps = function(chain, counts, group){
  chain@counts = as.integer(counts)
  chain@group = as.integer(group)

  chain@sums_n = as.integer(apply(counts, 2, sum))
  chain@sums_g = as.integer(apply(counts, 1, sum))

  chain@N = ncol(counts)
  chain@G = nrow(counts)

  if(!length(chain@samples_return)) chain@samples_return = sample.int(chain@N, 1)
  if(!length(chain@features_return)) chain@features_return = sample.int(chain@G, 1)

  if(!length(chain@samples_return_eps)) chain@samples_return_eps = sample.int(chain@N, 1)
  if(!length(chain@features_return_eps)) chain@features_return_eps = sample.int(chain@G, 1)

  for(n in paste0("samples_return", c("", "_eps")))
    slot(chain, n) = slot(chain, n)[1 <= slot(chain, n) & slot(chain, n) <= chain@N]

  for(n in paste0("features_return", c("", "_eps")))
    slot(chain, n) = slot(chain, n)[1 <= slot(chain, n) & slot(chain, n) <= chain@G]

  for(n in paste0(paste0(c("samples", "features"), "_return"), c("", "", "_eps", "_eps")))
    slot(chain, n) = sort(slot(chain, n))

  chain@Nreturn = length(chain@samples_return)
  chain@Greturn = length(chain@features_return)

  chain@NreturnEps = length(chain@samples_return_eps)
  chain@GreturnEps = length(chain@features_return_eps)

  for(n in c("hph", "lph", "mph"))
    slot(chain, n) = rep(0, chain@G)

  for(n in c(hyperparameters()))
    if(as.logical(chain@returns[n]))
      slot(chain, n) = rep(0, chain@M)

  if(as.logical(chain@returns["rho"])) chain@rho = rep(0, chain@M * chain@Nreturn)
  if(as.logical(chain@returns["eps"])) chain@eps = rep(0, chain@M * chain@NreturnEps * chain@GreturnEps)

  for(n in c("phi", "alp", "del", "gam", "xiPhi", "xiAlp", "xiDel"))
    if(as.logical(chain@returns[n]))
      slot(chain, n) = rep(0, chain@M * chain@Greturn)

  for(post in c("PostMean", "PostMeanSq")){
    for(s in paste0(hyperparameters(), post))
      slot(chain, s) = 0
    for(s in paste0(c("phi", "alp", "del", "gam", "xiPhi", "xiAlp", "xiDel"), post))
      slot(chain, s) = rep(0, chain@G)
    slot(chain, paste0("rho", post)) = rep(0, chain@N)
    slot(chain, paste0("eps", post)) = rep(0, chain@N * chain@G)
  }

  return(chain)
}
