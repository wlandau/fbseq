#' @title Function \code{scenario_heterosis_model}
#' @description produces a generic heterosis scenario
#' @export
#' @return a \code{Scenario} object for the heterosis problem
#' @param genes number of genes/genes in the data
#' @param libraries number of libraries/libraries in the data
#' @param design Design matrix, with \code{N} rows and \code{L} columns.
#' Same format as in \code{edgeR::glmFit}, with RNA-seq libraries/samples as rows
#' and each \code{beta_{l, .}} as a column.
#' @param truth a \code{Starts} object containing the necessary initialization constants
#' and hyperparameters.
scenario_heterosis_model = function(genes = 3.5e4, libraries = 12,
  design = cbind(rep(1, libraries), rep(c(1, -1, 1), each = floor(libraries/3)), rep(c(-1, 1, 1), each = floor(libraries/3))),
  truth = Starts(nu = 10, omegaSquared = 0.01, sigmaSquared = c(1, 0.25, 0.25), tau = 0.1, theta = c(3, 0, 0))){

  s = generate_data_from_model(genes = genes, design = design, truth = truth)

  s@bounds = rep(0, 4)
  names(s@bounds) = c("high-parent_1", "high-parent_2", "low-parent_1", "low-parent_2")

  s@contrasts = list(
    c(0, 1, 0),
    c(0, 0, 1),
    c(0, -1, 0),
    c(0, 0, -1))

  names(s@contrasts) = names(s@bounds)
  for(i in 1:length(s@contrasts)) names(s@contrasts[[i]]) = paste0("beta_", 1:ncol(s@design))

  rownames(s@design) = colnames(s@counts)
  colnames(s@design) = paste0("beta_", 1:ncol(s@design))

  s@propositions = list(1:2, 3:4)
  names(s@propositions) = c("high-parent_heterosis", "low-parent_heterosis")
  for(i in 1:length(s@propositions))
    names(s@propositions[[i]]) = names(s@contrasts)[s@propositions[[i]]]

  s
}
