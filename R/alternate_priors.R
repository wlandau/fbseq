#' @title Function \code{alternate_priors}
#' @details alternate priors for phi_g, alpha_g, and delta_g
#' @export
#' @return string vector of names of alternate priors
alternate_priors = function(){
  c("Laplace", "t", "horseshoe")
}