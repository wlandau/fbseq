#' @title Function \code{probs_single}
#' @description Extracts estimated posterior probs from a single \code{Chain} object. 
#'
#' @export
#' @return a data frame of estimated posterior probabilities, with rows corresponding to genes and columns 
#' corresponding to the propositions in \code{Scenario(chain)@@propositions}.
#' @param chain a \code{Chain} object
probs_single = function(chain){
  out = matrix(chain@probs, ncol = chain@P)
  rownames(out) = chain@gene_names
  colnames(out) = chain@proposition_names
  out
}

#' @title Function \code{probs}
#' @description Extracts estimated posterior probs from a \code{Chain} object.
#' or a list of \code{Chain} objects returned by \code{fbseq()}.
#'
#' @export
#' @return a data frame of estimated posterior probabilities, with rows corresponding to genes and columns 
#' corresponding to the propositions in the relevant \code{Scenario} object.
#' @param obj a \code{Chain} object or a list of \code{Chain} objects returned by \code{fbseq()}.
probs = function(obj){
  if(class(obj) == "Chain") obj = list(obj)
  pl = lapply(obj, function(ch) probs_single(ch) * ch@iterations * ch@thin)
  out = pl[[1]]
  if(length(pl) > 1) for(i in 2:length(pl)) out = out + pl[[i]]
  out/sum(sapply(obj, function(ch) ch@iterations * ch@thin))
}
