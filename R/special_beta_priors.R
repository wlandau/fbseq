#' @title Function \code{special_beta_priors}
#' @description special priors on the beta_{l, g} parameters (using the xi's as aux variables)
#' @export
#' @return string vector of names of special_beta priors
special_beta_priors = function(){
  c("Laplace", "t", "horseshoe")
}