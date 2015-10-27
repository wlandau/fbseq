#' @title Class \code{Starts}
#' @description A way for users to specify model parameter 
#' starting values. Use only when creating new \code{obj}
#' objects from scratch.
#' @exportClass Starts
#' 
#' @slot aRho initialization constant 
#' @slot aGamma initialization constant  
#' @slot bRho initialization constant 
#' @slot bGamma initialization constant 
#' @slot c initialization constants 
#' @slot dRho initialization constant 
#' @slot dGamma initialization constant 
#' @slot k initialization constants 
#' @slot r initialization constants 
#' @slot s initialization constants 
#' 
#' @slot beta MCMC starting values
#' @slot epsilon MCMC starting values
#' @slot gamma MCMC starting values
#' @slot nuGamma MCMC starting values
#' @slot nuRho MCMC starting values
#' @slot rho MCMC starting values
#' @slot sigmaSquared MCMC starting values
#' @slot tauGamma MCMC starting values
#' @slot tauRho MCMC starting values
#' @slot theta MCMC starting values
#' @slot xi MCMC starting values
setClass("Starts", 
  slots = list(
    aRho = "numeric",
    aGamma = "numeric",
    bRho = "numeric",
    bGamma = "numeric",
    c = "numeric",
    dRho = "numeric",
    dGamma = "numeric",
    k = "numeric",
    r = "numeric",
    s = "numeric",

    beta = "numeric",
    epsilon = "numeric",
    gamma = "numeric",
    nuGamma = "numeric",
    nuRho = "numeric",
    rho = "numeric",
    sigmaSquared = "numeric",
    tauGamma = "numeric",
    tauRho = "numeric",
    theta = "numeric",
    xi = "numeric"
  ),

  prototype = list(
    aRho = 1,
    aGamma = 1,
    bRho = 1,
    bGamma = 1,
    c = 10,
    dRho = 1000,
    dGamma = 1000,
    k = 5,
    r = 5,
    s = 100,

    tauGamma = 1
  )
)

#' @title Constructor for class \code{Starts}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' @export
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... additional slots.
Starts = function(obj = NULL, ...){
  starts = new("Starts", ...)

  if(class(obj) == "list") {
    for(n in slotNames(starts)){
      x = paste(n, "Start", sep = "")
      if(x %in% names(obj) && n %in% slotNames(starts))
        slot(starts, n) = obj[[x]]
      else if(n %in% intersect(names(obj), slotNames(starts)))
        slot(starts, n) = obj[[n]]
    }
  } else if(class(obj) == "Chain") {
    for(n in slotNames(starts)){
      x = paste(n, "Start", sep = "")
      if(x %in% slotNames(obj) && n %in% slotNames(starts))
        slot(starts, n) = slot(obj, x)
      else if(n %in% intersect(slotNames(obj), slotNames(starts)))
        slot(starts, n) = slot(obj, n)
    }
  }

  return(starts)
}