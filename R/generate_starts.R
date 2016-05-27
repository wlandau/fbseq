#' @title Function \code{get_nonzeros}
#' @description Helper function for \code{generate_starts()}.
#' @export
#' @return x with no zeros
#' @param x argument
get_nonzeros = function(x){
  if(!any(x > 0)) x = rexp(length(x))
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
  counts[counts < 1] = 0.5
  logcounts = log(counts)

  if(!length(starts@h)) starts@h = colMeans(logcounts) - mean(colMeans(logcounts))
  if(length(starts@h) != N) starts@h = rep(starts@h, length.out = N)
  h = matrix(rep(starts@h, each = G), ncol = N)
  logcounts = logcounts - h

  OLS = solve(t(design) %*% design) %*%  t(design)
  beta = t(OLS %*% t(logcounts))

  theta = apply(beta, 2, mean)
  sigmaSquared = apply(beta, 2, var)
  xi = rep(1, ncol(beta)*G)

  epsilon = logcounts - t(design %*% t(beta))
  gamma = get_nonzeros(apply(epsilon, 1, var, na.rm = T))

  nt = nu_tau(gamma)
  nu = nt$nu
  tau = nt$tau

  if(!is.finite(tau)) tau = 1e-12
  if(tau < 0) tau = 1e-12
  if(!is.finite(nu)) nu = 1e-12
  if(nu < 0) nu = 1e-12
  if(nu > starts@d) nu = starts@d

  for(n in c("c", "k", "q", "r", "s")){
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
#' @param upper upper bound for some parameters.
dispersed_set = function(chain, parm, lower = -Inf, upper = Inf){
  m = chain@iterations * chain@thin
  Mean = slot(chain, paste0(parm, "PostMean"))
  MeanSquare = slot(chain, paste0(parm, "PostMeanSquare"))
  Sd =  sqrt(m*(MeanSquare - Mean^2)/(m - 1))
  n = length(Mean)

  if(!length(lower)) lower = -Inf
  if(!length(upper)) upper = Inf

  lower = rep(lower, n/length(lower))
  upper = rep(upper, n/length(upper))

  lower = pmax(lower, Mean - 4*Sd)
  upper = pmin(upper, Mean + 4*Sd)

  df = 5
  out = rep(NA, n)
  ct = 0
  while(any(is.na(out))){
    ct = ct + 1
    if(ct > 1e3) stop("cannot disperse starting values. Try increasing iterations or thinning.")
    i = is.na(out)
    out[i] = rt(sum(i), df = df)*Sd[i]*sqrt((df-2)/df) + Mean[i]
    out[out <= lower & Sd > 0] = NA
    out[out >= upper & Sd > 0] = NA
  }
  out

#  Min = pmax(lower, Mean - 3*Sd, na.rm = T)
#  Max = pmin(upper, Mean + 3*Sd, na.rm = T) 
#  runif(n, Min, Max)
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
  chain@thin = as.integer(max(1, chain@thin))
  configs = Configs(chain)
  lower = list(nu = 0, gamma = min(Starts(chain)@gamma), sigmaSquared = 0, tau = 0, xi = min(Starts(chain)@xi))
  upper = list(nu = chain@d, sigmaSquared = chain@s^2)

  for(v in configs@parameter_sets_update)
    slot(chain, paste0(v, "Start")) = dispersed_set(chain, v, lower[[v]], upper[[v]])

  chain
}
