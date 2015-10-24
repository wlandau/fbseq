#' @title Function \code{concatenate}
#' @description Concatenate 2 \code{Chain} objects that are chronologically adjacent. 
#' @details If you run \code{chain1} and if \code{chain2 = run_mcmc(chain1)}, then you can call
#' \code{concatenate(chain1, chain2)} to combine them. The parameter samples 
#' will be combined, the starting values will be the same as for chain2, and 
#' \code{concatenate(chain1, chain2)@@M = chain1@@M + chain2@@M}. Also,
#' \code{concatenate(chain1, chain2)@@hph = (chain1@@M * chain1@@hph + chain2@@M * chain2@@hph) / (chain1@@M + chain2@@M)}
#' The other heterosis probabilities are similarly combined.
#' @export
#' @return a \code{Chain} object
#' @param chain1 a \code{Chain} object
#' @param chain2 a \code{Chain} object
concatenate = function(chain1, chain2){
  chain = chain2
  chain@M = chain1@M + chain2@M

  for(sn in "psrf")
    if(!length(slot(chain, sn)))
      slot(chain, sn) = slot(chain1, sn)

  sn = slotNames(chain)
  sn = sn[grepl("PostMean", sn)]
  sn = c(paste0(c("h", "l", "m"), "ph"), sn)

  for(p in sn)
    slot(chain, p) = (chain1@M * slot(chain1, p) + chain2@M * slot(chain2, p))/(chain1@M + chain2@M)

  for(p in parameters())
    if(chain@returns[p])
      slot(chain, p) = c(slot(chain1, p), slot(chain2, p))

  chain
}