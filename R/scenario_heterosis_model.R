#' @title Function \code{scenario_heterosis_model}
#' @description produces a generic heterosis scenario
#' @export
#' @return a \code{Scenario} object for the heterosis problem
#' @param genes number of genes/genes in the data
#' @param libraries number of libraries/libraries in the data
#' @param truth \code{Starts} object of hyperparameters underlying simulation
#' all hyperparameters must be specified.
scenario_heterosis_model = function(genes = 3e4, libraries = 16,
  truth = Starts(nu = 2.812224, 
    omegaSquared = 0.006808205, tau = 0.006780517, sigmaSquared = c(1, 
        0.03724313, 0.0324207, 0.0006229287, 0.06410533), theta = c(3, 
        -0.005734982, -0.02541216, -0.004763663, -0.06341044))){

  data(paschold)
  paschold = get("paschold")
  stopifnot(libraries >= nrow(paschold@design))

  ns = 0:(libraries -1) %% ncol(paschold@counts) + 1
  design = paschold@design[ns,]
  s = generate_data_from_model(genes = genes, design = design, truth = truth)

  libnames = colnames(paschold@counts)
  libnames = gsub("B73xMo17_Mo17xB73", "hybrids", libnames)
  libnames = gsub("B73xMo17", "hybrid1", libnames)
  libnames = gsub("Mo17xB73", "hybrid2", libnames)
  libnames = gsub("B73", "parent1", libnames)
  libnames = gsub("Mo17", "parent2", libnames)
  libnames = gsub("_.*", "", libnames)
  libnames = paste0(libnames[ns], "_", 1:libraries)

  colnames(s@counts) = libnames
  rownames(s@design) = libnames
  colnames(s@design) = paste0("beta_", 1:ncol(design))

  for(n in c("bounds", "contrasts", "propositions"))
    slot(s, n) = slot(paschold, n)

  cnames = names(s@bounds)
  cnames = gsub("B73xMo17_Mo17xB73", "hybrids", cnames)
  cnames = gsub("B73xMo17", "hybrid1", cnames)
  cnames = gsub("Mo17xB73", "hybrid2", cnames)
  names(s@bounds) = names(s@contrasts) = cnames
  for(i in 1:length(s@propositions))
    names(s@propositions[[i]]) = names(s@contrasts)[s@propositions[[i]]]

  pnames = names(s@propositions)
  pnames = gsub("B73xMo17_Mo17xB73", "hybrids", pnames)
  pnames = gsub("B73xMo17", "hybrid1", pnames)
  pnames = gsub("Mo17xB73", "hybrid2", pnames)
  names(s@propositions) = pnames

  s
}
