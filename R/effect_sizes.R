#' @title Function \code{effect_sizes}
#' @description Computes the "effect sizes" of propositions involving contrasts.
#'
#' @export
#' @return a data frame of estimated effect sizes, with rows corresponding to genes and columns 
#' corresponding to the propositions in \code{Scenario(chain)@@propositions}.
#' @param obj a \code{Chain} object or a list of \code{Chain} objects returned by \code{fbseq()}.
effect_sizes = function(obj){
  if(class(obj) == "Chain") obj = list(obj)
  est = estimates(obj)
  b = matrix(est[grep("beta_", rownames(est)), "mean"], nrow = obj[[1]]@G)
  s = Scenario(obj[[1]])
  out = NULL
  for(prop in s@propositions){
    eff = rep(Inf, obj[[1]]@G)
    for(c in prop){
      cst = b %*% s@contrasts[[c]] - s@bounds[c]
      cst = cst * (cst > 0)
      eff = pmin(eff, cst)
    }
    out = cbind(out, eff)
  }
  colnames(out) = obj[[1]]@proposition_names
  rownames(out) = obj[[1]]@gene_names
  out
}