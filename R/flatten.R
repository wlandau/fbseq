#' @title Function \code{flatten}
#' @description Flattens an object into some kind of useful array.
#' @export
#' 
#' @return one of the following, depending on the class of the argument object.
#' \describe{
#'   \item{\code{Chain} argument}{a data frame parameter libraries.}
#'   \item{\code{Starts} argument}{a named numeric vector of starting values. If
#' \code{epsilon} is set in \code{obj}, then one of \code{rho}, \code{phi}, \code{alp},
#' or \code{gamma} must be set too.}
#' }
#' @param obj The argument, either a \code{Chain} or a \code{Starts} object.
flatten = function(obj){
  if(class(obj) == "Chain")
    return(flatten_chain(obj))
  else if(class(obj) == "Starts")
    return(flatten_starts(obj))
  else
    stop("argument must be of class Chain or class Starts.")
}

#' @title Function \code{flatten_chain}
#' @description Flattens a \code{Chain} object.
#' @export
#' @return an array of MCMC parameter libraries.
#' @param chain a \code{Chain} object.
flatten_chain = function(chain){
  ret = list()

  for(x in c("nu", "omegaSquared", "tau"))
    if(length(slot(chain, x)))
      ret[[x]] = slot(chain, x)
  
  for(x in c("sigmaSquared", "theta"))
    if(length(slot(chain, x))){
      ret[[x]] = matrix(slot(chain, x), ncol = chain@L, byrow = T)
      colnames(ret[[x]]) = paste0(x, "_", 1:chain@L)
    }

  for(x in c("beta", "xi"))
    if(length(slot(chain, x))){
      ret[[x]] = matrix(slot(chain, x), ncol = chain@L * chain@Greturn, byrow = T)
      colnames(ret[[x]]) = paste0(x, "_", rep(1:chain@L, each = chain@Greturn), "_", rep(chain@genes_return, times = chain@L))
    }

  if(length(chain@epsilon)){
    ret$epsilon = matrix(chain@epsilon, ncol = chain@NreturnEpsilon * chain@GreturnEpsilon, byrow = T)
    colnames(ret$epsilon) = paste("epsilon_", rep(chain@libraries_return_epsilon, each = chain@GreturnEpsilon), 
      "_", rep(chain@genes_return_epsilon, times = chain@NreturnEpsilon), sep="")
  }

  if(length(chain@gamma)){
    ret$gamma = matrix(chain@gamma, ncol = chain@Greturn, byrow = T)
    colnames(ret$gamma) = paste("gamma", chain@genes_return, sep = "_")
  }

  if(length(chain@rho)){
    ret$rho = matrix(chain@rho, ncol = chain@Nreturn, byrow = T)
    colnames(ret$rho) = paste("rho", chain@libraries_return, sep = "_")
  }

  as.data.frame(do.call(cbind, ret))
}

#' @title Function \code{flatten_starts}
#' @description Flattens a \code{Starts} object.
#' @export
#' @return an array of MCMC parameter libraries.
#' @param starts a \code{Starts} object.
flatten_starts = function(starts){
  G = tryCatch(length(starts@gamma), warning = function(w) 0, error = function(w) 0)
  L = length(starts@beta)/ifelse(G, G, 1)
  N = tryCatch(max(length(starts@rho), length(starts@epsilon)/G), warning = function(w) 0, error = function(w) 0)

  for(x in c("beta", "xi"))
    if(length(slot(starts, x)))
      names(slot(starts, x)) = paste0(x, "_", rep(1:L, each = G), "_", rep(1:G, times = L))

  for(x in c("c", "k", "q", "r", "s", "gamma", "rho", "sigmaSquared", "theta"))
    if(length(slot(starts, x)))
      names(slot(starts, x)) = paste(x, 1:length(slot(starts, x)), sep = "_")

  if(length(starts@epsilon))
    names(starts@epsilon) = paste("epsilon_", rep(1:N, each = G), "_", 
                                                       rep(1:G, times = N), sep="")

  c(
    a = starts@a,
    b = starts@b,
    starts@c,
    d = starts@d,
    starts@k,
    starts@q,
    starts@r,
    starts@s,
    starts@w,

    starts@beta,
    starts@epsilon,
    starts@gamma,
    nu = starts@nu,
    omegaSquared = starts@omegaSquared,
    starts@rho,
    starts@sigmaSquared,
    tau = starts@tau,
    starts@theta,
    starts@xi
  )
}

#' @title Function \code{flatten_post}
#' @description Flattens posterior means or posterior means of squares
#' @export
#' @return an array of MCMC parameter libraries.
#' @param chain a \code{Chain} object.
#' @param square TRUE/FALSE value. If TRUE, return the flattened posterior mean squares
#' instead of posterior means
#' @param updated_only TRUE/FALSE. If TRUE, return the posterior quantities only for the 
#' parameters that were updated in the MCMC.
flatten_post = function(chain, square = F, updated_only = T){
  post = ifelse(square, "PostMeanSquare", "PostMean")
  faux_starts = Starts()
  for(x in slotNames(chain)[grep(paste0(post, "$"), slotNames(chain))])
    slot(faux_starts, gsub(post, "", x)) = slot(chain, x)

  u = chain@parameter_sets_update
  cand = parameters()
  if(updated_only) cand = intersect(cand, names(u)[as.logical(u)])
  pattern = paste0(cand, collapse = "|")
  flat = flatten(faux_starts)
  flat[grep(pattern, names(flat))]
}
