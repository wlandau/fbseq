#' @title Function \code{get_nonzeros}
#' @description Helper function for \code{generate_starts()}.
#' @export
#' @return x with no zeros
#' @param x argument
get_nonzeros = function(x){
  stopifnot(any(x > 0))
  x1 = x[x > 0]
  x[x <= 0] = sample(x1, sum(x <= 0), replace = T)
  x
}

#' @title Function \code{nu_tau}
#' @description Helper function for \code{generate_starts()}.
#' @export
#' @return list of nu and tau estimates
#' @param x argument
nu_tau = function(x){
  m = mean(x)
  v = var(x)
  shape = m^2/v + 2
  scale  = m*(shape - 1)
  nu = 2*shape
  tau = sqrt(2*scale/nu)
  list(nu = nu, tau = tau)
}

#' @title Function \code{generate_starts}
#' @description Generates MCMC starting values using a \code{Starts} object and an 
#' \code{Internals} object. The \code{configure(...)} function must be called first. 
#' @details Simple nonstochastic starting values are calculated by quick arithmetic. 
#' @export
#' 
#' @return a \code{Starts} object with a full set of MCMC starting values.
#' @param counts A data frame/matrix of RNA-seq read counts or list of slots. If a list,
#' then the function will return a \code{Chain} object with those slots.
#' @param design Gene-specific design matrix.
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
#' @param starts a \code{Starts} object to fill
generate_starts = function(counts, design, starts = Starts()){
  N = dim(counts)[2]
  G = dim(counts)[1]
  colnames(design) = NULL

  counts = as.matrix(counts)
  logcounts = log(counts + 1)

  OLS = solve(t(design) %*% design) %*%  t(design)
  PROJ = design %*% OLS
  beta = t(OLS %*% t(logcounts))

  theta = apply(beta, 2, mean)
  sigmaSquared = apply(beta, 2, var)
  xi = rep(1, ncol(beta)*G)

  epsilon = logcounts - t(PROJ %*% t(logcounts))
  rho = colMeans(epsilon)
  omegaSquared = var(rho)

  rhomat = matrix(rep(rho, each = G), ncol = N)
  gamma = get_nonzeros(apply(epsilon - rhomat, 1, var, na.rm = T))

  nt = nu_tau(gamma)
  nu = nt$nu
  tau = nt$tau

  for(n in c("c", "k", "r", "s")){
    if(length(slot(starts, n)) == 0)
      slot(starts, n) = slot(Starts(), n)
    if(length(slot(starts, n)) == 1)
      slot(starts, n) = rep(slot(starts, n), ncol(design))
    stopifnot(length(slot(starts, n)) == ncol(design))
  }

  for(n in slotNames(starts))
    if(!length(slot(starts, n)))
      slot(starts, n) = as(get(n), class(slot(starts, n)))

  starts
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
  if(!length(lower)) lower = NA
  if(!length(upper)) upper = NA
  m = chain@iterations
  Mean = slot(chain, paste0(parm, "PostMean"))
  MeanSquare = slot(chain, paste0(parm, "PostMeanSquare"))
  Sd =  sqrt(m*(MeanSquare - Mean^2)/(m - 1))
  n = length(Mean)

#  out = rt(n, 5)*Sd + Mean
#  if(all(is.finite(lower))) out[out <= lower] = lower
#  if(all(is.finite(upper))) out[out >= upper] = upper
#  out

  Min = pmax(rep(lower, n/length(lower)), Mean - 3*Sd, na.rm = T)
  Max = pmin(rep(upper, n/length(upper)), Mean + 3*Sd, na.rm = T) 
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
  lower = list(nu = 0, gamma = min(Starts(chain)@gamma), omegaSquared = 0, sigmaSquared = 0, 
                    tau = 0, xi = min(Starts(chain)@xi))
  upper = list(nu = chain@d, omegaSquared = chain@w^2, sigmaSquared = chain@s^2)

  for(v in configs@parameter_sets_update)
    slot(chain, paste0(v, "Start")) = dispersed_set(chain, v, lower[[v]], upper[[v]])

  chain
}
