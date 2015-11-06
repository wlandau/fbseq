#' @title Function \code{scenario_heterosis_model}
#' @description produces a generic heterosis scenario
#' @export
#' @return a \code{Scenario} object for the heterosis problem
#' @param genes number of genes/genes in the data
#' @param libraries number of libraries/libraries in the data
scenario_heterosis_model = function(genes = 3.5e4, libraries = 16){
  stopifnot(!(libraries %% 8))

  truth = Starts(nu = 10, omegaSquared = 0.01, tau = 0.1,
                        sigmaSquared = c(1, 0.25, 0.25, 0.25, 0.25), 
                        theta = c(3, 0, 0, 0, 0))

  design = cbind(
    rep(1, libraries/4), 
    rep(c(1, -1, 1, 1), each = libraries/4), 
    rep(c(-1, 1, 1, 1), each = libraries/4),
    rep(c(0, 0, 1, -1), each = libraries/4),
    rep(rep(c(1, -1), each = ceiling(libraries/8)), times = 4))

  colnames(design) = paste0("beta_", 1:5)
  rownames(design) = paste0("library", 1:libraries)

  s = generate_data_from_model(genes = genes, design = design, truth = truth)

  s@contrasts = list(
    "high-parent_hybrids_1" = c(beta_1 = 0, beta_2 =  1, beta_3 =  0, beta_4 = 0, beta_5 = 0),
    "high-parent_hybrids_2" = c(beta_1 = 0, beta_2 =  0, beta_3 =  1, beta_4 = 0, beta_5 = 0),
    "low-parent_hybrids_1"  = c(beta_1 = 0, beta_2 = -1, beta_3 =  0, beta_4 = 0, beta_5 = 0),
    "low-parent_hybrids_2"  = c(beta_1 = 0, beta_2 = 0,  beta_3 = -1, beta_4 = 0, beta_5 = 0),

    "high-parent_hybrid1_1" = c(beta_1 = 0, beta_2 =  0, beta_3 =  2, beta_4 =  1, beta_5 = 0),
    "high-parent_hybrid1_2" = c(beta_1 = 0, beta_2 =  2, beta_3 =  0, beta_4 =  1, beta_5 = 0),
    "low-parent_hybrid1_1"  = c(beta_1 = 0, beta_2 =  0, beta_3 = -2, beta_4 = -1, beta_5 = 0),
    "low-parent_hybrid1_2"  = c(beta_1 = 0, beta_2 = -2, beta_3 =  0, beta_4 = -1, beta_5 = 0),

    "high-parent_hybrid2_1" = c(beta_1 = 0, beta_2 =  0, beta_3 =  2, beta_4 = -1, beta_5 = 0),
    "high-parent_hybrid2_2" = c(beta_1 = 0, beta_2 =  2, beta_3 =  0, beta_4 = -1, beta_5 = 0),
    "low-parent_hybrid2_1"  = c(beta_1 = 0, beta_2 =  0, beta_3 = -2, beta_4 =  1, beta_5 = 0),
    "low-parent_hybrid2_2"  = c(beta_1 = 0, beta_2 = -2, beta_3 =  0, beta_4 =  1, beta_5 = 0))

  s@bounds = rep(0, length(s@contrasts))
  names(s@bounds) = names(s@contrasts)

  s@propositions = list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12)
  names(s@propositions) = gsub("_2", "", names(s@contrasts)[1:6*2])

  for(i in 1:length(s@propositions))
    names(s@propositions[[i]]) = names(s@contrasts)[s@propositions[[i]]]

  s
}
