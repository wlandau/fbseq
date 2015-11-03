#' @title Class \code{Inference}
#' @description A collection of MCMC control parameters.
#' @exportClass Inference
#' 
#' @slot conjunctions List of numeric vectors of length J. Each numeric vector is a collection of 
#' indices c_1, ..., c_J such that the posterior probability,
#' P(contrasts[[c_1]] * beta > value[c_1], ..., contrasts[[c_J]] * beta > value[c_J] | data),
#' is estimated by averging over Monte Carlo samples. This allows the user to "test"
#' multiple contrasts together. See the vignette for a more thorough explanation.
#' @slot contrasts list of numeric vectors, each length \code{L}. Each vector
#' is a contrast of the \code{beta_{l, .}}'s of interest.
#' @slot design Design matrix, with \code{N} rows and \code{L} columns.
#' Same format as in \code{edgeR::glmFit}, with RNA-seq libraries/samples as rows
#' and each \code{beta_{l, .}} as a column.
#' @slot probs matrix with G rows and J columns. These are the posterior probabilities
#' of each conjunction above for each gene, estimated by averaging over Monte Carlo
#' samples.
#' @slot values Numeric vector of length \code{length(contrasts)}. Used int tests.
setClass("Inference", 
  slots = list(
    conjunctions = "list",
    contrasts = "list",
    design = "matrix",
    probs = "matrix",
    values = "numeric"
  )
)

#' @title Constructor for class \code{Inference}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' @export
#' @return a \code{Inference} object
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... optional slot values.
Inference = function(obj = NULL, ...){
  inference = new("Inference", ...)

  if(class(obj) == "list") {
    for(n in intersect(slotNames(inference), names(obj)))
       slot(inference, n) = obj[[n]]
  } else if(class(obj) == "Chain"){
    stop("not implemented yet.")
  }

  return(Inference)
}