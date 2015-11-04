#' @name paschold
#' @title paschold
#' @description a \code{Scenario} object guiding the analysis of the Paschold dataset.
#' @docType data
#' @author Paschold
NULL

#' @name starting_chain
#' @title starting_chain
#' @description \code{Chain} object for use in the vignette.
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
NULL

#' @name ending_chain
#' @title ending_chain
#' @description \code{Chain} object for use in the vignette.
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
NULL

#' @name example_generated_data
#' @title example_generated_data
#' @description Example generated data generated with \code{generate_data()}.
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
NULL

#' @title Function \code{make_vignette_data}
#' @description Creates the vignette data.
#' @export
make_vignette_data = function(){
  file = "vignette.rda"
  example_generated_data = dat = generate_data(libraries = 12, genes = 20)
  starting_chain = Chain(dat$counts, dat$design, Configs(diag = "none", ess = 0, iterations = 100, burnin = 100, thin = 0))
  ending_chain = fbseq(starting_chain)
  save(example_generated_data, starting_chain, ending_chain, file = file)
}
