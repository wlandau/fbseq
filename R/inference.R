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
    design = "matrix",
    contrasts = "list",
    values = "numeric",
    conjunctions = "list",
    probs = "matrix"
  )
)

#' @title Function \code{check_inference}
#' @details checks an \code{Inference} object for inconsistencies.
#' @export
#' @return a \code{Inference} object
#' @param inference a \code{Inference} object
check_inference = function(inference){
  for(s in slotNames(inference)) stopifnot(all(is.finite(unlist(slot(inference, s)))))
  if(any(!dim(inference@design))) stop("no design matrix specified.")

  N = nrow(inference@design)
  L = ncol(inference@design)

  if(!all(sapply(inference@contrasts, length) == L)) stop("the number of terms in each contrast must equal the number of columns in the design matrix.")
  if(length(inference@contrasts) != length(inference@values)) 
    stop("the length of the \"values\" vector must equal the number of contrasts.")
  if(!all(unique(unlist(inference@conjunctions)) %in% 1:length(inference@contrasts))) 
    stop("elements of the \"conjunctions\" list must be vectors whose elements lie beween 1 and the number of contrasts.")
  if(length(inference@conjunctions) & !length(inference@contrasts)) 
    stop("\"conjunctions\" list specified without a list of contrasts.")
  if(length(inference@conjunctions) & !length(inference@values)) 
    stop("\"conjunctions\" list specified without a \"values\" vector.")  

  inference
}

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

  check_inference(inference)
}