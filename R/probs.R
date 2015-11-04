#' @title Function \code{probs}
#' @description Extracts estimated posterior probs from a \code{Chain} object. 
#'
#' @export
#' @return a data frame of estimated posterior probabilities, with rows corresponding to genes and columns 
#' corresponding to the propositions in \code{Scenario(chain)@@propositions}.
#' @param chain a \code{Chain} object
probs = function(chain){
  out = matrix(chain@probs, ncol = chain@P)
  rownames(out) = chain@gene_names
  colnames(out) = chain@proposition_names
  out
}
