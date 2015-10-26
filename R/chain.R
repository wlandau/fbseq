#' @title Class \code{Chain}
#' @description the main storage object for a heterosis MCMC chain.
#' Create a new one with a call to \code{Chain()} and run MCMC
#' by feeding one into \code{single_chain()}.
#' @exportClass Chain
#'
#' @slot diag convergence diagnostic to use. Can be "gelman" or "none".
#' @slot ess Minimum effective sample size for all parameters
#' @slot max_attempts Maximum number of retries for assessing convergence and generating enough effective samples.
#' Can be set to Inf to run indefinitely.
#' @slot nchains_diag number of independent chains to run (including this one) to use
#' convergence diagnostics that require multiple chains. 
#' @slot psrf_tol upper threshold for Gelman-Rubin potential scale reduction factors (if diag is "gelman")
#' 
#' @slot burnin MCMC burnin, the number of MCMC iterations to ignore at the beginning of each obj
#' @slot genes_return Indices of genes whose parameter samples you want to return.
#' Applies to all gene-specific parameters except for the epsilons.
#' @slot genes_return_epsilon Indices of genes g for which epsilon_{n, g} is updated/returned.
#' @slot iterations Number of MCMC iterations (not including burnin or thinning)
#' @slot libraries_return Indices of RNA-seq libraries whose parameter samples you want to return.
#' @slot libraries_return_epsilon Indices of RNA-seq libraries n for which epsilon_{n, g} is updated/returned.
#' Applies to all library-specific parameters except for the epsilons.
#' @slot parameter_sets_return Character vector naming the variables whose MCMC samples 
#' you want to return
#' @slot parameter_sets_update Character vector naming the variables to calculate/update
#' during the MCMC.
#' @slot priors Names of the family of priors on the betas after integrating out the xi's. 
#' Can be any value returned by \code{alternate_priors()}. All other values will default to the normal prior.
#' @slot thin MCMC thinning interval, number of iterations to skip in between iterations to return.
#' @slot verbose Number of times to print out progress during burnin and the actual MCMC.
#' If \code{verbose} > 0, then progress messages will also print during setup and cleanup.
#' 
#' @slot counts RNA-seq count data, flattened from a matrix
#' @slot countSums_g gene-specific count sums
#' @slot countSums_n library-specific count sums
#' @slot design Gene-specific design, flattened from the design matrix. Must contain only 0's, 1's, and -1's.
#' Original matrix must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
#' @slot G number of genes
#' @slot Greturn number of genes to return gene-specific MCMC parameter samples for (except the epsilons)
#' @slot GreturnEpsilon number of genes to return gene-specific MCMC epsilon parameter samples
#' @slot L number of columns in the original design matrix
#' @slot N number of libraries
#' @slot Nreturn number of libraries to return library-specific MCMC parameter samples for (except the epsilons)
#' @slot NreturnEpsilon number of libraries to return library-specific MCMC epsilon parameter samples
#' @slot psrf Gelman-Rubin potential scale reduction factors for the sampled parameters (even if the 
#' actual MCMC parameter samples are not returned)
#' @slot seeds vector of N*G random number generator seeds
#' 
#' @slot aRho initialization constant 
#' @slot aGamma initialization constant  
#' @slot bRho initialization constant 
#' @slot bGamma initialization constant 
#' @slot c initialization constants 
#' @slot dRho initialization constant 
#' @slot dGamma initialization constant 
#' @slot k initialization constants 
#' @slot r initialization constants 
#' @slot s initialization constants 
#' 
#' @slot beta MCMC parameter samples
#' @slot epsilon MCMC parameter samples
#' @slot gamma MCMC parameter samples
#' @slot nuGamma MCMC parameter samples
#' @slot nuRho MCMC parameter samples
#' @slot omega MCMC parameter samples
#' @slot rho MCMC parameter samples
#' @slot tauGamma MCMC parameter samples
#' @slot tauRho MCMC parameter samples
#' @slot theta MCMC parameter samples
#' @slot xi MCMC parameter samples
#' 
#' @slot betaStart MCMC starting values
#' @slot epsilonStart MCMC starting values
#' @slot gammaStart MCMC starting values
#' @slot nuRhoStart MCMC starting values
#' @slot nuGammaStart MCMC starting values
#' @slot omegaStart MCMC starting values
#' @slot rhoStart MCMC starting values
#' @slot tauRhoStart MCMC starting values
#' @slot tauGammaStart MCMC starting values
#' @slot thetaStart MCMC starting values
#' @slot xiStart MCMC starting values
#' 
#' @slot betaPostMean estimated posterior means
#' @slot epsilonPostMean estimated posterior means
#' @slot gammaPostMean estimated posterior means
#' @slot nuRhoPostMean estimated posterior mean
#' @slot nuGammaPostMean estimated posterior mean
#' @slot omegaPostMean estimated posterior means
#' @slot rhoPostMean estimated posterior means
#' @slot tauRhoPostMean estimated posterior mean
#' @slot tauGammaPostMean estimated posterior mean
#' @slot thetaPostMean estimated posterior means
#' @slot xiPostMean estimated posterior means
#' 
#' @slot betaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot epsilonPostMeanSquare estimated posterior means of the squares of parameters
#' @slot gammaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot nuRhoPostMeanSquare estimated posterior mean of the square of the parameter
#' @slot nuGammaPostMeanSquare estimated posterior mean of the square of the parameter 
#' @slot omegaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot rhoPostMeanSquare estimated posterior means of the squares of parameters
#' @slot tauRhoPostMeanSquare estimated posterior mean of the square of the parameter
#' @slot tauGammaPostMeanSquare estimated posterior mean of the square of the parameter
#' @slot thetaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot xiPostMeanSquare posterior means of the squares of parameters
setClass("Chain",
  slots = list(
    diag = "character",
    ess = "integer",
    max_attempts = "numeric",
    nchains_diag = "integer",
    psrf_tol = "numeric",

    burnin = "integer",
    genes_return = "integer",
    genes_return_epsilon = "integer",
    iterations = "integer",
    libraries_return = "integer",
    libraries_return_epsilon = "integer",
    parameter_sets_return = "integer",
    parameter_sets_update = "integer",
    priors = "integer",
    thin = "integer",
    verbose = "integer",

    counts = "integer",
    countSums_g = "integer",
    countSums_n = "integer",
    design = "integer",
    G = "integer",
    Greturn = "integer",
    GreturnEpsilon = "integer",
    L = "integer",
    N = "integer",
    Nreturn = "integer",
    NreturnEpsilon = "integer",
    psrf = "numeric",
    seeds = "integer",

    aRho = "numeric",
    aGamma = "numeric",
    bRho = "numeric",
    bGamma = "numeric",
    c = "numeric",
    dRho = "numeric",
    dGamma = "numeric",
    k = "numeric",
    r = "numeric",
    s = "numeric",

    beta = "numeric",
    epsilon = "numeric",
    gamma = "numeric",
    nuGamma = "numeric",
    nuRho = "numeric",
    omega = "numeric",
    rho = "numeric",
    tauGamma = "numeric",
    tauRho = "numeric",
    theta = "numeric",
    xi = "numeric",

    betaStart = "numeric",
    epsilonStart = "numeric",
    gammaStart = "numeric",
    nuRhoStart = "numeric",
    nuGammaStart = "numeric",
    omegaStart = "numeric",
    rhoStart = "numeric",
    tauRhoStart = "numeric",
    tauGammaStart = "numeric",
    thetaStart = "numeric",
    xiStart = "numeric",

    betaPostMean = "numeric",
    epsilonPostMean = "numeric",
    gammaPostMean = "numeric",
    nuRhoPostMean = "numeric",
    nuGammaPostMean = "numeric",
    omegaPostMean = "numeric",
    rhoPostMean = "numeric",
    tauRhoPostMean = "numeric",
    tauGammaPostMean = "numeric",
    thetaPostMean = "numeric",
    xiPostMean = "numeric",

    betaPostMeanSquare = "numeric",
    epsilonPostMeanSquare = "numeric",
    gammaPostMeanSquare = "numeric",
    nuRhoPostMeanSquare = "numeric",
    nuGammaPostMeanSquare = "numeric",
    omegaPostMeanSquare = "numeric",
    rhoPostMeanSquare = "numeric",
    tauRhoPostMeanSquare = "numeric",
    tauGammaPostMeanSquare = "numeric",
    thetaPostMeanSquare = "numeric",
    xiPostMeanSquare = "numeric"
  )
)
