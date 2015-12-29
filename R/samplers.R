#' @include sampler_options.R
NULL

#' @title Class \code{Samplers}
#' @description specification of the MCMC black box sampler of each parameter
#' @exportClass Samplers
#' 
#' @slot beta sampler specification (character vector of length 1, element of samplers())
#' @slot epsilon sampler specification (character vector of length 1, element of samplers())
#' @slot gamma sampler specification (character vector of length 1, element of samplers())
#' @slot nu sampler specification (character vector of length 1, element of samplers())
#' @slot sigmaSquared sampler specification (character vector of length 1, element of samplers())
#' @slot tau sampler specification (character vector of length 1, element of samplers())
#' @slot theta sampler specification (character vector of length 1, element of samplers())
#' @slot xi sampler specification (character vector of length 1, element of samplers())
setClass("Samplers",
  slots = list(
    beta = "character",
    epsilon = "character",
    gamma = "character",
    nu = "character",
    sigmaSquared = "character",
    tau = "character",
    theta = "character",
    xi = "character"
  ),
  prototype = list(
    beta = "default",
    epsilon = "default",
    gamma = "default",
    nu = "default",
    sigmaSquared = "default",
    tau = "default",
    theta = "default",
    xi = "default"
  )
)

#' @title Constructor for class \code{Samplers}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' @export
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... additional slots.
Samplers = function(obj = NULL, ...){
  samplers = new("Samplers", ...)
  if(length(obj) == 1 && class(obj) == "character")
    for(n in slotNames(samplers))
      slot(samplers, n) = obj

  if(class(obj) == "list") {
    for(n in slotNames(samplers)){
      x = paste(n, "Sampler", sep = "")
      if(x %in% names(obj) && n %in% slotNames(samplers))
        slot(samplers, n) = as(obj[[x]], class(slot(samplers, n)))
      else if(n %in% intersect(names(obj), slotNames(samplers)))
        slot(samplers, n) = as(obj[[n]], class(slot(samplers, n)))
    }
  } else if(class(obj) == "Chain") {
    for(n in slotNames(samplers)){
      x = paste(n, "Sampler", sep = "")
      slot(samplers, n) = sampler_options()[slot(obj, x)]
    }
  }

  for(n in slotNames(samplers)) stopifnot(slot(samplers, n) %in% sampler_options() && length(slot(samplers, n)) == 1)
  samplers
}
