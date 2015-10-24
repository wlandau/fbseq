#' @name paschold_counts
#' @title paschold_counts
#' @description Paschold count data matrix.
#' @docType data
#' @author Paschold
NULL

#' @name paschold_group
#' @title paschold_group
#' @description A vector of integers denoting the experimental design.
#' There is one integer for each RNA-seq sample/library, each denoting
#'  the genetic variety of that sample: 1 for parent 1, 2 for the hybrid,
#' and 3 for parent 2. 
#' @docType data
#' @author Will Landau \email{will.landau@@gmail.com}
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
  file = "vignette_data.rda"
  example_generated_data = dat = generate_data(samples = 12, features = 20)
  starting_chain = Chain(dat$counts, dat$group, Configs(diag = "none", ess = 0, max_attempts = 5, iterations = 100, burnin = 100, thin = 0))
  ending_chain = heterosis(starting_chain)
  save(example_generated_data, starting_chain, ending_chain, file = file)
}
