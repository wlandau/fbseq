#' @title Function \code{smallchain}
#' @description Make a small \code{Chain} object for debugging and testing.
#' @export
#' @return A small \code{Chain} object
#' @param size Number of genes. Also the number of iterations and burnin.
smallchain = function(size = 25){
  data(paschold, envir = environment())
  paschold = get("paschold")
  size = 25
  paschold@counts = paschold@counts[1:size,]
  con = Configs(iterations = size, burnin = size, thin = 1, priors = "Laplace")
  Chain(paschold, con)
}
