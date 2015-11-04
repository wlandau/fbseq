#' @title Function \code{effect_sizes}
#' @description Computes the "effect sizes" of propositions involving contrasts.
#'
#' @export
#' @return a data frame of estimated effect sizes, with rows corresponding to genes and columns 
#' corresponding to the propositions in \code{Scenario(chain)@@propositions}.
#' @param chain a \code{Chain} object
effect_sizes = function(chain){
  b = matrix(chain@betaPostMean, nrow = chain@G)
  s = Scenario(chain)
  out = NULL
  for(prop in s@propositions){
    eff = rep(Inf, chain@G)
    for(c in prop){
      cst = b %*% s@contrasts[[c]] - s@bounds[c]
      cst = cst * (cst > 0)
      eff = pmin(eff, cst)
    }
    out = cbind(out, eff)
  }
  colnames(out) = chain@proposition_names
  rownames(out) = chain@gene_names
  out
}