#' @title Function \code{flatten}
#' @description Flattens an object into some kind of useful array.
#' @export
#' 
#' @return one of the following, depending on the class of the argument object.
#' \describe{
#'   \item{\code{Chain} argument}{a data frame parameter samples.}
#'   \item{\code{Starts} argument}{a named numeric vector of starting values. If
#' \code{eps} is set in \code{obj}, then one of \code{rho}, \code{phi}, \code{alp},
#' or \code{gam} must be set too.}
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
#' @return an array of MCMC parameter samples.
#' @param chain a \code{Chain} object.
flatten_chain = function(chain){
  ret = list()

  for(x in hyperparameters())
    if(length(slot(chain, x)))
      ret[[x]] = slot(chain, x)
  
  if(length(chain@rho)){
    ret$rho = matrix(chain@rho, ncol = chain@Nreturn, byrow = T)
    colnames(ret$rho) = paste("rho", chain@samples_return, sep = "_")
  }
  
  for(x in c("phi", "alp", "del", "gam", "xiPhi", "xiAlp", "xiDel"))
    if(length(slot(chain, x))){
      ret[[x]] = matrix(slot(chain, x), ncol = chain@Greturn, byrow = T)
      colnames(ret[[x]]) = paste(x, chain@features_return, sep = "_")
    }
  
  if(length(chain@eps)){
    ret$eps = matrix(chain@eps, ncol = chain@NreturnEps * chain@GreturnEps, byrow = T)
    colnames(ret$eps) = paste("eps_", rep(chain@samples_return_eps, each = chain@GreturnEps), 
      "_", rep(chain@features_return_eps, times = chain@NreturnEps), sep="")
  }

  as.data.frame(do.call(cbind, ret))
}

#' @title Function \code{flatten_starts}
#' @description Flattens a \code{Starts} object.
#' @export
#' @return an array of MCMC parameter samples.
#' @param starts a \code{Starts} object.
flatten_starts = function(starts){

  G = tryCatch(max(length(starts@phi), length(starts@alp), length(starts@del), length(starts@gam)),
    warning = function(w) 0, error = function(w){0})
  N = tryCatch(max(length(starts@rho), length(starts@eps)/G), warning = function(w) 0, error = function(w){0})

  if(!is.finite(N)) N = 0
  if(!is.finite(G)) G = 0

  if(N && !G) G = tryCatch(length(starts@eps)/N, warning = function(w) 0, error = function(w){0})
  if(length(starts@eps) && !(N || G))
    stop("cannot call flatten(starts) if eps is set but not rho, phi, alp, del, gam, xiPhi, xiAlp, or xiDel.")

  for(x in c("phi", "alp", "del", "rho", "gam", "xiPhi", "xiAlp", "xiDel"))
    if(length(slot(starts, x)))
      names(slot(starts, x)) = paste(x, 1:length(slot(starts, x)), sep = "_")

  if(length(starts@eps))
    names(starts@eps) = paste("eps_", rep(1:N, each = G), "_", 
                                                       rep(1:G, times = N), sep="")

  c(
    dRho = starts@dRho,
    dGam = starts@dGam,
    aRho = starts@aRho,
    aGam = starts@aGam,
    bRho = starts@bRho,
    bGam = starts@bGam,
    cPhi = starts@cPhi,
    cAlp = starts@cAlp,
    cDel = starts@cDel,
    sPhi = starts@sPhi,
    sAlp = starts@sAlp,
    sDel = starts@sDel,
    kPhi = starts@kPhi,
    kAlp = starts@kAlp,
    kDel = starts@kDel,
    rPhi = starts@rPhi,
    rAlp = starts@rAlp,
    rDel = starts@rDel,
    nuRho = starts@nuRho,
    nuGam = starts@nuGam,
    tauRho = starts@tauRho,
    tauGam = starts@tauGam,
    thePhi = starts@thePhi,
    theAlp = starts@theAlp,
    theDel = starts@theDel,
    sigPhi = starts@sigPhi,
    sigAlp = starts@sigAlp,
    sigDel = starts@sigDel,
    starts@phi,
    starts@alp,
    starts@del,
    starts@rho,
    starts@gam,
    starts@xiPhi,
    starts@xiAlp,
    starts@xiDel,
    starts@eps
  )
}

#' @title Function \code{flatten_post}
#' @description Flattens posterior means or posterior means of squares
#' @export
#' @return an array of MCMC parameter samples.
#' @param chain a \code{Chain} object.
#' @param square TRUE/FALSE value. If TRUE, return the flattened posterior mean squares
#' instead of posterior means
#' @param updated_only TRUE/FALSE. If TRUE, return the posterior quantities only for the 
#' parameters that were updated in the MCMC.
flatten_post = function(chain, square = F, updated_only = T){
  post = ifelse(square, "PostMeanSq", "PostMean")
  N = chain@N
  G = chain@G

  for(v in hyperparameters())
    names(slot(chain, paste0(v, post))) = v

  vs = c("phi", "alp", "del", "gam", "xiPhi", "xiAlp", "xiDel")
  for(v in vs)
    names(slot(chain, paste0(v, post))) = paste0(v, "_", 1:G)

  names(slot(chain, paste0("rho", post))) = paste0("rho_", 1:N)
  names(slot(chain, paste0("eps", post))) = paste("eps_", rep(1:N, each = G), "_", rep(1:G, times = N), sep="")

  out = NULL
  cand = c(hyperparameters(), "rho", vs, "eps")
  if(updated_only) cand = setdiff(cand, names(which(!chain@updates)))
  for(v in cand) out = c(out, slot(chain, paste0(v, post)))

  out
}
