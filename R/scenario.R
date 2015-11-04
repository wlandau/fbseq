#' @title Class \code{Scenario}
#' @description A collection of MCMC control parameters.
#' @exportClass Scenario
#' 
#' @slot conjunctions List of numeric vectors of length J. Each numeric vector is a collection of 
#' indices c_1, ..., c_J such that the posterior probability,
#' P(contrasts[[c_1]] * beta > value[c_1], ..., contrasts[[c_J]] * beta > value[c_J] | data),
#' is estimated by averging over Monte Carlo samples. This allows the user to "test"
#' multiple contrasts together. See the vignette for a more thorough explanation.
#' @slot contrasts list of numeric vectors, each length \code{L}. Each vector
#' is a contrast of the \code{beta_{l, .}}'s of interest.
#' @slot counts RNA-seq count data matrix with G rows (for genes) and N columns (for libraries/samples).
#' @slot design Design matrix, with \code{N} rows and \code{L} columns.
#' Same format as in \code{edgeR::glmFit}, with RNA-seq libraries/samples as rows
#' and each \code{beta_{l, .}} as a column.
#' @slot values Numeric vector of length \code{length(contrasts)}. Used int tests.
setClass("Scenario", 
  slots = list(
    conjunctions = "list",
    contrasts = "list",
    counts = "matrix",
    design = "matrix",
    values = "numeric"
  )
)

#' @title Function \code{check_scenario}
#' @details checks an \code{Scenario} object for inconsistencies.
#' @export
#' @return a \code{Scenario} object
#' @param scenario a \code{Scenario} object
check_scenario = function(scenario){
  for(s in slotNames(scenario)) stopifnot(all(is.finite(unlist(slot(scenario, s)))))

  if(any(!dim(scenario@counts))) stop("no count data specified.")
  if(any(!dim(scenario@design))) stop("no design matrix specified.")
  if(ncol(scenario@counts) != nrow(scenario@design)) stop("ncol(counts) must equal nrow(design).")

  designUniqueN = as.integer(apply(scenario@design, 2, function(x){length(unique(x[x != 0]))}))
  if(any(!designUniqueN)) stop("every column in the design matrix must have nonzero elements.")

  N = nrow(scenario@design)
  L = ncol(scenario@design)

  if(!all(sapply(scenario@contrasts, length) == L)) stop("the number of terms in each contrast must equal the number of columns in the design matrix.")
  if(length(scenario@contrasts) != length(scenario@values)) 
    stop("the length of the \"values\" vector must equal the number of contrasts.")
  if(!all(unique(unlist(scenario@conjunctions)) %in% 1:length(scenario@contrasts))) 
    stop("elements of the \"conjunctions\" list must be vectors whose elements lie beween 1 and the number of contrasts.")
  if(length(scenario@conjunctions) & !length(scenario@contrasts)) 
    stop("\"conjunctions\" list specified without a list of contrasts.")
  if(length(scenario@conjunctions) & !length(scenario@values)) 
    stop("\"conjunctions\" list specified without a \"values\" vector.")  

  scenario
}

#' @title Constructor for class \code{Scenario}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' Elements passed with \code{...} must be named. For example, \code{Scenario(design = my_matrix)}.
#' @export
#' @return a \code{Scenario} object
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... optional slot values.
Scenario = function(obj = NULL, ...){
  scenario = new("Scenario", ...)

  if(class(obj) == "list") {
    for(n in intersect(slotNames(scenario), names(obj)))
       slot(scenario, n) = obj[[n]]
  } else if(class(obj) == "Chain"){
    cm = matrix(obj@conjunctions, nrow = obj@J)
    scenario@conjunctions = lapply(1:dim(cm)[2], function(i) which(cm[,i] == 1))
    scenario@contrasts = as.list(as.data.frame(matrix(obj@contrasts, nrow = obj@L)))
    scenario@counts = matrix(obj@counts, ncol = obj@N)
    scenario@design = matrix(obj@design, nrow = obj@N)
    scenario@values = obj@values
  }

  check_scenario(scenario)
}