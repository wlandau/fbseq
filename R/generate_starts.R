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

  if(tau < 0) tau = 1e-12
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
  n = chain@iterations * chain@thin
  Mean = flatten_post(chain, square = F, updated_only = T)
  MeanSq = flatten_post(chain, square = T, updated_only = T)
  Sd =  sqrt(n*(MeanSq - Mean^2)/(n - 1))

  iter = 100
  con0 = Configs(chain)
  con = Configs(
    iterations = iter,
    burnin = con0@thin,
    thin = con0@thin,
    priors = con0@priors,
    parameter_sets_return = con0@parameter_sets_update,
    parameter_sets_update = con0@parameter_sets_update,
    genes_return = 1:chain@G,
    genes_return_epsilon = 1,
    libraries_return = 1:chain@N,
    libraries_return_epsilon = 1,
    verbose = con0@verbose)

    if(chain@verbose) print("Running a mini chain to disperse starting values.")
    mini = Chain(Scenario(chain), con)
    mini = single_mcmc(mini)
    samples = mcmc_samples(mini)
    ns = intersect(names(Mean), colnames(samples))
    ns = ns[!grepl("epsilon_", ns)]
    samples = samples[,ns]
    Mean = Mean[ns]
    Sd = Sd[ns]
    
    scaled = sweep(samples, 2, Mean)
    scaled = sweep(scaled, 2, Sd, "/")
    norms = rowSums(scaled^2)
    i = sample((1:iter)[order(norms)][ceiling(0.95*iter):iter], 1)
    starts = as.numeric(samples[i,])
    names(starts) = colnames(samples)

    out = chain
    parms = names(starts)
    for(k in 1:10) parms = gsub("_[0-9]*", "", parms)
    for(n in intersect(con0@parameter_sets_update, parms))
      slot(out, paste0(n, "Start")) = starts[n == parms]
    out
}
