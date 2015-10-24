#' @title Function \code{nu_tau}
#' @description Helper function for \code{simple_starts()}.
#' @export
#' @return list of nu and tau estimates
#' @param x argument
nu_tau = function(x){
  m = mean(x^2)
  v = var(x^2)
  shape = m^2/v + 2
  scale  = m*(shape - 1)
  nu = 2*shape
  tau = sqrt(2*scale/nu)
  list(nu = nu, tau = tau)
}

#' @title Function \code{simple_starts}
#' @description Generates MCMC starting values using a \code{Starts} object and an 
#' \code{Internals} object. The \code{configure(...)} function must be called first. 
#' @details Simple nonstochastic starting values are calculated by quick arithmetic. 
#' @export
#' 
#' @return a \code{Chain} object with a full set of MCMC starting values.
#' @param chain \code{Chain} object whose starting values to fill.
#' @param counts A data frame/matrix of RNA-seq read counts or list of slots. If a list,
#' then the function will return a \code{Chain} object with those slots.
#' @param group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for parent 2,
#' and 3 for the hybrid.
simple_starts = function(chain, counts, group){
  N = dim(counts)[2]
  G = dim(counts)[1]

  counts = as.matrix(counts)
  logcounts = log(counts + 1)

  rotation = 0.5 * cbind(
    c(1, 1, 0), 
    c(0, -1, 1), 
    c(-1, 0, 1))

  n = table(group)[as.character(1:3)]
  mu = sapply(1:3, function(i) rowMeans(matrix(logcounts[, group == i], ncol = n[i])))
  rotated = mu %*% rotation

  phi = rotated[,1]
  alp = rotated[,2]
  del = rotated[,3]

  thePhi = mean(phi)
  theAlp = mean(alp)
  theDel = mean(del)

  sigPhi = sd(phi)
  sigAlp = sd(alp)
  sigDel = sd(del)

  xiPhi = rep(1, G)
  xiAlp = rep(1, G)
  xiDel = rep(1, G)

  eps = logcounts - mu[, group]

  gam = apply(eps, 1, sd, na.rm = T)
  if(all(gam == 0)) gam = 1/rgamma(G, 1)
  if(any(gam == 0)) gam[gam == 0] = min(gam[gam != 0])
  gammat = matrix(rep(gam, times = N), ncol = N)

  rho = apply(eps/gammat, 2, sd, na.rm = T)
  if(all(rho == 0)) rho = 1/rgamma(N, 10)
  if(any(rho == 0)) rho[rho == 0] = min(rho[rho != 0])

  eps = as.numeric(eps)

  nt = nu_tau(rho)
  nuRho = nt$nu
  tauRho = nt$tau

  nt = nu_tau(gam)
  nuGam = nt$nu
  tauGam = nt$tau

  for(n in slotNames(chain))
    if(grepl("Start", n) && !length(slot(chain, n)))
      slot(chain, n) = get(gsub("Start", "", n))

  return(chain)
}


#' @title Function \code{dispersed_set}
#' @description Produces a dispersed starting values (relative to the estimated posterior so far)
#' for a given set of parameters.
#' 
#' @export
#' @rdname dispersed_set
#' @aliases dispersed_set
#' @return a vector with a dispersed starting values (relative to the estimated posterior so far)
#' for a given set of parameters.
#' @param chain \code{Chain} object that has already been run with \code{run_mcmc()}.
#' @param parm character string specifying the parameter set to work on.
#' @param lower lower bound for nonnegative parameters.
#' @param upper upper bound for nonnegative parameters.
dispersed_set = function(chain, parm, lower = NA, upper = NA){
  m = chain@M
  Mean = slot(chain, paste0(parm, "PostMean"))
  MeanSq = slot(chain, paste0(parm, "PostMeanSq"))
  Sd =  sqrt(m*(MeanSq - Mean^2)/(m - 1))
  n = length(Mean)
  Min = pmax(rep(lower, n), Mean - 3*Sd, na.rm = T)
  Max = pmin(rep(upper, n), Mean + 3*Sd, na.rm = T) 
#  ifelse(runif(n) < 0.5, Min, Max)
  runif(n, Min, Max)
}

#' @title Function \code{disperse_starts}
#' @description Uses the hyperparameter samples from a previous chain to 
#' generate another \code{Chain} object with dispersed starting values.
#' Multiple \code{Chain} objects produced by \code{disperse_startes()} can
#' be run and Gelman potential scale reduction factors on the collective output
#' will make sense as convergence diagnostics
#' 
#' @export
#' @rdname disperse_starts
#' @aliases disperse_starts
#' @return a \code{Chain} object with dispersed MCMC starting values.
#' @param chain \code{Chain} object that has already been run with \code{run_mcmc()}.
disperse_starts = function(chain){
  configs = Configs(chain)
  lower = c(nuRho = 0, nuGam = 0, tauRho = 0, tauGam = 0, sigPhi = 0, sigAlp = 0, sigDel = 0,
                  rho = 0, gam = 0, xiPhi = 0, xiAlp = 0, xiDel = 0)
  upper = c(nuRho = chain@dRho, nuGam = chain@dGam, 
                  sigPhi = chain@sPhi, sigAlp = chain@sAlp, sigDel = chain@sDel)

  for(v in configs@updates)
    slot(chain, paste0(v, "Start")) = dispersed_set(chain, v, lower[v], upper[v])

  chain
}
