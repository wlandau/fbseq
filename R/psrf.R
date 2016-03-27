#' @include flatten.R
NULL

#' @title Function \code{psrf}
#' @description calculate all Gelman potential scale reduction factors on a list of parallel independent chains.
#' @export
#' @return Gelman potential scale reduction factors on a list of parallel independent chains.
#' @param chain_list list of \code{Chain} objects
psrf = function(chain_list){
  stopifnot(length(chain_list) > 1)
  stopifnot(length(unique(sapply(chain_list, function(chain) chain@iterations))) == 1)

  Mean = sapply(chain_list, flatten_post)
  MeanSquare = sapply(chain_list, flatten_post, square = T)

  m = length(chain_list)
  n = chain_list[[1]]@iterations

  SjSq =  n*(MeanSquare - Mean^2)/(n - 1)
  W = rowMeans(SjSq)
  B = n*apply(Mean, 1, var)
  sqrt((n-1 + B/W)/n)
}
