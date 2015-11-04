#' @include scenario_heterosis_model.R
NULL

#' @name paschold
#' @title paschold
#' @description a \code{Scenario} object guiding the analysis of the Paschold dataset.
#' @docType data
#' @author Paschold
NULL

#' @name chain1
#' @title chain1
#' @description \code{Chain} object for use in the vignette.
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
NULL

#' @name chain2
#' @title chain2
#' @description \code{Chain} object for use in the vignette.
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
NULL

#' @name tiny
#' @title tiny
#' @description Example \code{Scenario} object from \code{scenario_heterosis_model}
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
NULL

#' @title Function \code{make_vignette_data}
#' @description Creates the vignette data.
#' @export
make_vignette_data = function(){
  tiny = scenario_heterosis_model(genes = 20)
  chain1 = Chain(tiny, Configs(diag = "none", ess = 0, iterations = 100, burnin = 100, thin = 0))
  chain2 = fbseq(chain1)
  save(tiny, chain1, chain2, file = "vignette.rda")
}
