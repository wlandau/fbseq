#' @title Class \code{Scenario}
#' @description A collection of MCMC control parameters.
#' @exportClass Scenario
#' 
#' @slot bounds Numeric vector of length \code{length(contrasts)}. Used to formulate 
#' logical propositions involving contrasts.
#' @slot contrasts list of numeric vectors, each length \code{L}. Each vector
#' is a contrast of the \code{beta_{l, .}}'s of interest.
#' @slot counts RNA-seq count data matrix with G rows (for genes) and N columns (for libraries/samples).
#' @slot design Design matrix, with \code{N} rows and \code{L} columns.
#' Same format as in \code{edgeR::glmFit}, with RNA-seq libraries/samples as rows
#' and each \code{beta_{l, .}} as a column.
#' @slot propositions List of numeric vectors of length \code{P}. Each numeric vector \code{v} denotes a logical
#' proposition that involves the \code{contrasts} list and \code{bounds} vector. Each element of \code{v}
#' is the index of a contrast in \code{contrasts}. For each gene, the MCMC will estimate the posterior
#' probability that each contrast (of the \code{beta} parameters) denoted in \code{v} is greater than 
#' its corresponding value in \code{bound}. See the tutorial vignette for more details
setClass("Scenario", 
  slots = list(
    bounds = "numeric",
    contrasts = "list",
    counts = "matrix",
    design = "matrix",
    propositions = "list"
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

  L = ncol(scenario@design)
  N = nrow(scenario@design)

  if(!all(sapply(scenario@contrasts, length) == L)) stop("the number of terms in each contrast must equal the number of columns in the design matrix.")
  if(length(scenario@contrasts) != length(scenario@bounds))
    stop("the length of the \"bounds\" vector must equal the number of contrasts.")
  if(!all(unique(unlist(scenario@propositions)) %in% 1:length(scenario@contrasts)))
    stop("elements of the \"propositions\" list must be vectors whose elements lie beween 1 and the number of contrasts.")
  if(length(scenario@propositions) & !length(scenario@contrasts))
    stop("\"propositions\" list specified without a list of contrasts.")
  if(length(scenario@propositions) & !length(scenario@bounds))
    stop("\"propositions\" list specified without a \"bounds\" vector.")

  scenario
}

#' @title Constructor for class \code{Scenario}
#' @details Precedence will be given to the \code{Chain} or \code{list} object over \code{...}.
#' Elements passed with \code{...} must be named. For example, \code{Scenario(design = my_matrix)}.
#' @export
#' @return a \code{Scenario} object
#' @param obj a \code{Chain} or \code{list} object to get slots from.
#' @param ... optional slot bounds.
Scenario = function(obj = NULL, ...){
  scenario = new("Scenario", ...)

  if(class(obj) == "list") {
    for(n in intersect(slotNames(scenario), names(obj)))
       slot(scenario, n) = as(obj[[n]], class(slot(scenario, n)))
  } else if(class(obj) == "Chain"){
    scenario@bounds = obj@bounds
    scenario@contrasts = as.list(as.data.frame(matrix(obj@contrasts, nrow = obj@L)))
    scenario@counts = matrix(obj@counts, ncol = obj@N)
    scenario@design = matrix(obj@design, nrow = obj@N)
    cm = matrix(obj@propositions, ncol = obj@P)
    scenario@propositions = lapply(1:dim(cm)[2], function(i) which(cm[,i] == 1))

    names(scenario@bounds) = obj@bound_names
    names(scenario@contrasts) = obj@contrast_names
    colnames(scenario@counts) = obj@library_names
    rownames(scenario@counts) = obj@gene_names
    colnames(scenario@design) = paste0("beta_", 1:ncol(scenario@design))
    rownames(scenario@design) = obj@library_names
    names(scenario@propositions) = obj@proposition_names

    for(i in 1:length(scenario@contrasts))
      names(scenario@contrasts[[i]]) = paste0("beta_", 1:ncol(scenario@design))

    for(i in 1:length(scenario@propositions))
      names(scenario@propositions[[i]]) = names(scenario@contrasts)[scenario@propositions[[i]]]
  }

  check_scenario(scenario)
}