#' @title Function \code{parameters}
#' @description Prints out the full collection of variable names. See the methodology vignette for details.
#'
#' @export
#' @return Character vector of parameter names.
parameters = function(){
  c("beta", "epsilon", "gamma", "nu", "omegaSquared", "rho", "sigmaSquared", "tau", "theta", "xi")
}

#' @title Function \code{contstants}
#' @description Prints out initialization constants of the model. See the methodology vignette for details.
#'
#' @export
#' @return Character vector of parameter names.
constants = function(){
  c("a", "b", "c", "d", "k", "r", "s", "w")
}
