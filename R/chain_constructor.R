#' @include alternate_priors.R generate_starts.R
NULL

#' @title Function \code{Chain}
#' @description Construct a \code{Chain} object, the 
#' one argument to \code{single_chain(...)}.
#' @details If \code{list} is set, all other arguments will be dropped.
#' @export
#' @return a \code{Chain} object ready to be passed to \code{single_chain()}
#'
#' @param data A data frame or matrix of RNA-seq read counts or list of slots. If a list,
#' then the function will return a \code{Chain} object with those slots. Note: \code{data}
#' should only be a list within internal functions of the package. It is not recommended that
#' the user assign \code{data} to be a list.
#' @param inference an \code{Inference} object.
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting values.
Chain = function(
  data, 
  inference, 
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
  chain = plug_in_chain(chain, inference, configs = Configs(), starts = Starts())
  chain = plug_in_chain(chain, inference, configs = configs, starts = starts)
  chain = fill_easy_gaps(chain, counts, inference)
  chain = simple_starts(chain, counts, inference@design)
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
#' @param inference an \code{Inference} object
#' @param configs A \code{Configs} object of MCMC control parameters.
#' @param starts A \code{Starts} object of model parameter starting values.
plug_in_chain = function(chain, inference, configs, starts){
  chain@iterations = as.integer(configs@iterations)
  subtract = c("parameter_sets_return", "parameter_sets_update", "priors")
  for(n in setdiff(intersect(slotNames(chain), slotNames(configs)), subtract))
    slot(chain, n) = as(slot(configs, n), class(slot(chain, n)))

  for(n in c("parameter_sets_return", "parameter_sets_update")){
    slot(chain, n) = as.integer(parameters() %in% slot(configs, n))
    names(slot(chain, n)) = parameters()
  }

  chain@priors = ifelse(configs@priors %in% alternate_priors(), which(alternate_priors() == configs@priors), as.integer(0))
  if(length(chain@priors) == 1) chain@priors = rep(chain@priors, ncol(inference@design))
  stopifnot(length(chain@priors) == ncol(inference@design))

  for(n in slotNames(chain)){
    if(grepl("Start", n))
      slot(chain, n) = as(slot(starts, gsub("Start", "", n)), class(slot(chain, n)))
    else if(n %in% slotNames(starts) && !(paste0(n, "Start") %in% slotNames(chain)))
      slot(chain, n) = as(slot(starts, n), class(slot(chain, n)))
  }

  chain@conjunctions = as.integer(unlist(lapply(inference@conjunctions, function(x){1:length(inference@contrasts) %in% x})))
  chain@contrasts = unlist(inference@contrasts)
  chain@design = as.numeric(inference@design)
  chain@values = inference@values

  chain
}

#' @title Function \code{fill_easy_gaps}
#' @description Set required slots in a \code{Chain} object.
#' Used internally.
#' @export
#' @return a \code{Chain} object 
#'
#' @param chain a \code{Chain} object
#' @param counts Matrix of RNA-seq read counts.
#' @param inference an \code{Inference} object
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' sets of gene-specific variables.
fill_easy_gaps = function(chain, counts, inference){
  stopifnot(ncol(counts) == nrow(inference@design))

  designUnique = apply(inference@design, 2, function(x){
    out = sort(unique(x[x != 0]))
    c(out, rep(0, dim(inference@design)[1] - length(out)))
  })

  designUniqueN = as.integer(apply(inference@design, 2, function(x){length(unique(x[x != 0]))}))
  if(any(!designUniqueN)) stop("every column in the inference matrix must have nonzero elements.")

  if(!length(chain@betas_update)) chain@betas_update = 1:ncol(inference@design)

  chain@C = length(inference@contrasts)
  chain@counts = as.integer(counts)
  chain@countSums_g = as.integer(apply(counts, 1, sum))
  chain@countSums_n = as.integer(apply(counts, 2, sum))
  chain@designUnique = as.numeric(designUnique)
  chain@designUniqueN = as.integer(designUniqueN)
  chain@genes_return = sort(chain@genes_return)
  chain@genes_return_epsilon = sort(chain@genes_return_epsilon)
  chain@G = G = nrow(counts)
  chain@Greturn = Greturn = length(chain@genes_return)
  chain@GreturnEpsilon = GreturnEpsilon = length(chain@genes_return_epsilon)
  chain@libraries_return = sort(chain@libraries_return)
  chain@libraries_return_epsilon = sort(chain@libraries_return_epsilon)
  chain@J = length(inference@conjunctions)
  chain@L = L = ncol(inference@design)
  chain@Lupdate = length(chain@betas_update)
  chain@N = N = nrow(inference@design)
  chain@Nreturn = Nreturn = length(chain@libraries_return)
  chain@NreturnEpsilon = NreturnEpsilon = length(chain@libraries_return_epsilon)
  chain@probs = rep(0, chain@J * chain@G)

  lengths = c(
    beta = L*Greturn,
    epsilon = NreturnEpsilon*GreturnEpsilon,
    gamma = Greturn,
    nu = 1,
    omegaSquared = 1,
    rho = Nreturn,
    sigmaSquared = L,
    tau = 1,
    theta = L,
    xi = L*Greturn)

  for(s in names(lengths))
    if(as.logical(chain@parameter_sets_return[s]))
      slot(chain, s) = rep(0, chain@iterations*lengths[s])

  lengths = c(
    beta = L*G,
    epsilon = N*G,
    gamma = G,
    nu = 1,
    omegaSquared = 1,
    rho = N,
    sigmaSquared = L,
    tau = 1,
    theta = L,
    xi = L*G)

  for(post in c("PostMean", "PostMeanSquare"))
    for(s in names(lengths)){
      n = paste0(s, post)
      slot(chain, n) = rep(0, lengths[s])
    }

  chain
}
