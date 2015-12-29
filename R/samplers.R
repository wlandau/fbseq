#' @title Function \code{samplers}
#' @description Returns the vector of options for slot of the \code{Samplers} object. See the tutorial vignette for details.
#' @export
#' @return vector of options for slot of the \code{Samplers} object
samplers = function(){
  c("default", "slice_step", "random_walk_metropolis")
}
