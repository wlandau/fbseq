#' @title Class \code{Starts}
#' @description A way for users to specify model parameter 
#' starting values. Use only when creating new \code{obj}
#' objects from scratch.
#' @exportClass Starts
#' 
#' @slot dRho initialization constant
#' @slot dGam initialization constant
#' @slot aRho initialization constant
#' @slot aGam initialization constant
#' @slot bRho initialization constant
#' @slot bGam initialization constant
#' @slot cPhi initialization constant
#' @slot cAlp initialization constant
#' @slot cDel initialization constant
#' @slot sPhi initialization constant
#' @slot sAlp initialization constant
#' @slot sDel initialization constant
#' @slot kPhi initialization constant
#' @slot kAlp initialization constant
#' @slot kDel initialization constant
#' @slot rPhi initialization constant
#' @slot rAlp initialization constant
#' @slot rDel initialization constant
#' @slot nuRho model parameter 
#' @slot nuGam model parameter 
#' @slot tauRho model parameter
#' @slot tauGam model parameter 
#' @slot thePhi model parameter 
#' @slot theAlp model parameter 
#' @slot theDel model parameter 
#' @slot sigPhi model parameter 
#' @slot sigAlp model parameter 
#' @slot sigDel model parameter 
#' @slot phi model parameter 
#' @slot alp model parameter 
#' @slot del model parameter
#' @slot rho model parameter
#' @slot gam model parameter 
#' @slot xiPhi model parameter
#' @slot xiAlp model parameter
#' @slot xiDel model parameter
#' @slot eps model parameter 
setClass("Starts", 
  slots = list(
    dRho = "numeric",
    dGam = "numeric",
    aRho = "numeric",
    aGam = "numeric", 
    bRho = "numeric",
    bGam = "numeric",
    cPhi = "numeric",
    cAlp = "numeric",
    cDel = "numeric",
    sPhi = "numeric",
    sAlp = "numeric",
    sDel = "numeric",
    kPhi = "numeric",
    kAlp = "numeric",
    kDel = "numeric",
    rPhi = "numeric",
    rAlp = "numeric",
    rDel = "numeric",

    nuRho = "numeric",
    nuGam = "numeric",
    tauRho = "numeric",
    tauGam = "numeric",

    thePhi = "numeric",
    theDel = "numeric",
    theAlp = "numeric",

    sigPhi = "numeric",
    sigAlp = "numeric",
    sigDel = "numeric", 

    phi = "numeric",
    alp = "numeric",
    del = "numeric",
    rho = "numeric", 
    gam = "numeric",
    xiPhi = "numeric", 
    xiAlp = "numeric", 
    xiDel = "numeric", 
    eps = "numeric"
  ),

  prototype = list(
    dRho = 1e3,
    dGam = 1e3,
    aRho = 2, 
    aGam = 2, 
    bRho = 1,
    bGam = 1,
    cPhi = 10,
    cAlp = 10,
    cDel = 10,
    sPhi = 1000,
    sAlp = 1000,
    sDel = 1000,
    kPhi = 5,
    kAlp = 5,
    kDel = 5,
    rPhi = 5,
    rAlp = 5,
    rDel = 5,
    tauGam = 1
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