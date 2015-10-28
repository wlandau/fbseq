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
#' @param design Gene-specific design matrix.
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting values.
Chain = function(
  data, 
  design, 
  configs = Configs(),
  starts = Starts()
){
  chain = new("Chain")

  if(class(data) == "list" & !is.data.frame(data)){
    for(n in intersect(names(data), slotNames(chain)))
       slot(chain, n) = data[[n]]
    return(chain)
  }

  counts = as.matrix(data)
  stopifnot(all(is.finite(counts)))
  chain = plug_in_chain(chain, design, configs = Configs(), starts = Starts())
  chain = plug_in_chain(chain, design, configs = configs, starts = starts)
  chain = fill_easy_gaps(chain, counts, design)
  chain = simple_starts(chain, counts, design)
  if(chain@parameter_sets_return["tauGamma"] || chain@parameter_sets_update["tauGamma"]) warning("TauGamma should be set constant at 1. Sampling tauGamma will prevent the MCMC from converging, and returning even a constant tauGamma will muddle effective sample size diagnostics. Control tauGamma with the returns, updates, returns_skip, and updates_skip slots in your Configs object. See the package vignettes for details.")
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
#' @param design Gene-specific design matrix. 
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting values.
plug_in_chain = function(chain, design, configs, starts){
  chain@iterations = as.integer(configs@iterations)
  subtract = c("parameter_sets_return", "parameter_sets_update", "priors")
  for(n in setdiff(intersect(slotNames(chain), slotNames(configs)), subtract))
    slot(chain, n) = as(slot(configs, n), class(slot(chain, n)))

  for(n in c("parameter_sets_return", "parameter_sets_update")){
    slot(chain, n) = as.integer(parameters() %in% slot(configs, n))
    names(slot(chain, n)) = parameters()
  }

  chain@priors = ifelse(configs@priors %in% alternate_priors(), which(alternate_priors() == configs@priors), as.integer(0))
  if(length(chain@priors) == 1) chain@priors = rep(chain@priors, ncol(design))
  stopifnot(length(chain@priors) == ncol(design))

  for(n in slotNames(chain)){
    if(grepl("Start", n))
      slot(chain, n) = as(slot(starts, gsub("Start", "", n)), class(slot(chain, n)))
    else if(n %in% slotNames(starts) && !(paste0(n, "Start") %in% slotNames(chain)))
      slot(chain, n) = as(slot(starts, n), class(slot(chain, n)))
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
#' @param design Gene-specific design matrix. 
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
fill_easy_gaps = function(chain, counts, design){
  stopifnot(ncol(counts) == nrow(design))

  designUnique = apply(design, 2, function(x){
    out = sort(unique(x[x != 0]))
    c(out, rep(0, dim(design)[1] - length(out)))
  })

  designUniqueN = as.integer(apply(design, 2, function(x){length(unique(x[x != 0]))}))
  if(any(!designUniqueN)) stop("every column in the design matrix must have nonzero elements.")

  chain@counts = as.integer(counts)
  chain@countSums_g = as.integer(apply(counts, 1, sum))
  chain@countSums_n = as.integer(apply(counts, 2, sum))
  chain@design = as.numeric(design)
  chain@designUnique = as.numeric(designUnique)
  chain@designUniqueN = as.integer(designUniqueN)
  chain@G = G = nrow(counts)
  chain@Greturn = Greturn = length(chain@genes_return)
  chain@GreturnEpsilon = GreturnEpsilon = length(chain@genes_return_epsilon)
  chain@L = L = ncol(design)
  chain@N = N = nrow(design)
  chain@Nreturn = Nreturn = length(chain@libraries_return)
  chain@NreturnEpsilon = NreturnEpsilon = length(chain@libraries_return_epsilon)

  lengths = c(
    beta = L*Greturn,
    epsilon = NreturnEpsilon*GreturnEpsilon,
    gamma = Greturn,
    nuGamma = 1,
    nuRho = 1,
    rho = Nreturn,
    sigmaSquared = L,
    tauGamma = 1,
    tauRho = 1,
    theta = L,
    xi = L*Greturn)

  for(s in names(lengths))
    if(as.logical(chain@parameter_sets_return[s]))
      slot(chain, s) = rep(0, chain@iterations*lengths[s])

  lengths = c(
    beta = L*G,
    epsilon = N*G,
    gamma = G,
    nuGamma = 1,
    nuRho = 1,
    rho = N,
    sigmaSquared = L,
    tauGamma = 1,
    tauRho = 1,
    theta = L,
    xi = L*G)

  for(post in c("PostMean", "PostMeanSquare"))
    for(s in names(lengths)){
      n = paste0(s, post)
      slot(chain, n) = rep(0, lengths[s])
    }

  chain
}
