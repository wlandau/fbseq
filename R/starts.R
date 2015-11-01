#' @title Class \code{Starts}
#' @description A way for users to specify model parameter 
#' starting values. Use only when creating new \code{obj}
#' objects from scratch.
#' @exportClass Starts
#' 
#' @slot a initialization constant  
#' @slot b initialization constant 
#' @slot c initialization constants 
#' @slot d initialization constant 
#' @slot h log-scale normalization factors
#' @slot k initialization constants 
#' @slot r initialization constants 
#' @slot s initialization constants 
#' @slot w initialization constants 
#' 
#' @slot beta MCMC starting values
#' @slot delta MCMC starting values
#' @slot epsilon MCMC starting values
#' @slot gamma MCMC starting values
#' @slot nu MCMC starting values
#' @slot omegaSquared MCMC starting values
#' @slot pi MCMC starting values
#' @slot rho MCMC starting values
#' @slot sigmaSquared MCMC starting values
#' @slot tau MCMC starting values
#' @slot theta MCMC starting values
#' @slot xi MCMC starting values
setClass("Starts", 
  slots = list(
    a = "numeric",
    b = "numeric",
    c = "numeric",
    d = "numeric",
    h = "numeric",
    k = "numeric",
    r = "numeric",
    s = "numeric",
    w = "numeric",

    beta = "numeric",
    delta = "numeric",
    epsilon = "numeric",
    gamma = "numeric",
    nu = "numeric",
    omegaSquared = "numeric",
    pi = "numeric",
    rho = "numeric",
    sigmaSquared = "numeric",
    tau = "numeric",
    theta = "numeric",
    xi = "numeric"
  ),

  prototype = list(
    a = 1,
    b = 1,
    c = 10,
    d = 1000,
    k = 5,
    r = 5,
    s = 100,
    w = 100
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