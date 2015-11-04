#' @title Function \code{probs}
#' @description Extracts estimated posterior probs from a \code{Chain} object. 
#'
#' @export
#' @return a data frame of estimated posterior probabilities, with rows corresponding to genes and columns 
#' corresponding to the conjunctions in \code{Scenario(chain)@@conjunctions}.
#' @param chain a \code{Chain} object
probs = function(chain){
  out = matrix(chain@probs, ncol = chain@J)
  rownames(probs) = chain@gene_names
  colnames(probs) = chain@proposition_names
  out
}
