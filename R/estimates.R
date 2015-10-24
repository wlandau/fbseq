#' @include credible_intervals.R
NULL

#' @title Function \code{estimates}
#' @description Extracts heterosis probabilities, posterior means, and estimated credible intervals
#' from a \code{Chain} object. The order of the rows respects the order of the rows/features/genes
#' in the original RNA-seq count data.
#'
#' @export
#' @return A list of (1) "probs" a data frame of heterosis probabilities (with rows matching those of the 
#' original count data), (2) "params" a data frame of estimates: posterior means and approximate credible intervals.
#' Parameters not updeted in the MCMC are excluded, and (3) "psrf" Gelman-Rubin potential scale reduction factors
#' (empty if not using Gelman diagnostics)
#' @param chain a \code{Chain} object
#' @param level level of the credible intervals from 0 to 1
estimates = function(chain, level = 0.95){
  list(probs = data.frame(hph = chain@hph, lph = chain@lph, mph = chain@mph),
        params = credible_intervals(chain, level),
        psrf = chain@psrf)
}
