#' @title Function \code{assign_samplers}
#' @description Assigns the right samplers to the "Sampler" slots of a \code{Chain} object
#' @export
#' @return vector \code{Chain} object with samplers assigned
#' @param chain \code{Chain} object
#' @param configs \code{Configs} object
assign_samplers = function(chain, configs){
  s = configs@samplers
  ns = paste0(parameters(), "Sampler")

  if(s == "slice_step") {
    for(n in ns) slot(chain, n) = which(samplers() == "slice_step")
    chain@thetaSampler = which(samplers() == "default")
  } else if(s == "metropolis"){
    for(n in ns) slot(chain, n) = which(samplers() == "metropolis")
    chain@thetaSampler = which(samplers() == "default")    
  } else {
    for(n in ns) slot(chain, n) = which(samplers() == "default")
  }

  return(chain)
}