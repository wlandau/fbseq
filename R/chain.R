#' @title Class \code{Chain}
#' @description the main storage object for a heterosis MCMC chain.
#' Create a new one with a call to \code{Chain()} and run MCMC
#' by feeding one into \code{single_chain()}.
#' @exportClass Chain
#'
#' @slot psrf Gelman-Rubin potential scale reduction factors for the sampled parameters (even if the 
#' actual MCMC parameter samples are not returned)
#'
#' @slot conjunctions conjunctions of inequalities involving contrasts from the \code{Scenario} object.
#' @slot contrasts contrasts from the \code{Scenario} object.
#' @slot counts RNA-seq count data, flattened from a matrix
#' @slot design Design, flattened from the design matrix.
#' Original matrix must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' sets of gene-specific variables.
#' @slot values values to compare contrasts to (from the \code{Scenario} object)
#' 
#' @slot diag convergence diagnostic to use. Can be "gelman" or "none".
#' @slot ess Minimum effective sample size for all parameters
#' @slot max_attempts Maximum number of retries for assessing convergence and generating enough effective samples.
#' Can be set to Inf to run indefinitely.
#' @slot nchains_diag number of independent chains to run (including this one) to use
#' convergence diagnostics that require multiple chains. 
#' @slot psrf_tol upper threshold for Gelman-Rubin potential scale reduction factors (if diag is "gelman")
#' 
#' @slot betas_update values of l for which to update the beta_{l, g} parameters. Manually set for debugging purposes only.
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
#' @slot C number of contrasts
#' @slot countSums_g gene-specific count sums
#' @slot countSums_n library-specific count sums
#' @slot designUnique Matrix of unique nonzero elements of \code{design}. Vacent entries are 0.
#' @slot designUniqueN for each column index \code{l}, number of unique nonzero elements of \code{design[, l]}.
#' @slot G number of genes
#' @slot Greturn number of genes to return gene-specific MCMC parameter samples for (except the epsilons)
#' @slot GreturnEpsilon number of genes to return gene-specific MCMC epsilon parameter samples
#' @slot J number of conjunctions of contrasts
#' @slot L number of columns in the original design matrix
#' @slot Lupdate number of values of l for which to update the beta_{l, g} parameters. Manually set for debugging purposes only.
#' @slot N number of libraries
#' @slot Nreturn number of libraries to return library-specific MCMC parameter samples for (except the epsilons)
#' @slot NreturnEpsilon number of libraries to return library-specific MCMC epsilon parameter samples
#' @slot probs estimated posterior probabilities of conjunctions in the \code{Scenario} object.
#' @slot seeds vector of N*G random number generator seeds
#' 
#' @slot a initialization constant  
#' @slot b initialization constant 
#' @slot c initialization constants 
#' @slot d initialization constant 
#' @slot k initialization constants 
#' @slot r initialization constants 
#' @slot s initialization constants 
#' @slot w initialization constants 
#' 
#' @slot beta MCMC parameter samples
#' @slot epsilon MCMC parameter samples
#' @slot gamma MCMC parameter samples
#' @slot nu MCMC parameter samples
#' @slot omegaSquared MCMC parameter samples
#' @slot rho MCMC parameter samples
#' @slot sigmaSquared MCMC parameter samples
#' @slot tau MCMC parameter samples
#' @slot theta MCMC parameter samples
#' @slot xi MCMC parameter samples
#' 
#' @slot betaStart MCMC starting values
#' @slot epsilonStart MCMC starting values
#' @slot gammaStart MCMC starting values
#' @slot nuStart MCMC starting values
#' @slot omegaSquaredStart MCMC starting values
#' @slot rhoStart MCMC starting values
#' @slot sigmaSquaredStart MCMC starting values
#' @slot tauStart MCMC starting values
#' @slot thetaStart MCMC starting values
#' @slot xiStart MCMC starting values
#' 
#' @slot betaPostMean estimated posterior means
#' @slot epsilonPostMean estimated posterior means
#' @slot gammaPostMean estimated posterior means
#' @slot nuPostMean estimated posterior mean
#' @slot omegaSquaredPostMean estimated posterior means
#' @slot rhoPostMean estimated posterior means
#' @slot sigmaSquaredPostMean estimated posterior means
#' @slot tauPostMean estimated posterior mean
#' @slot thetaPostMean estimated posterior means
#' @slot xiPostMean estimated posterior means
#' 
#' @slot betaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot epsilonPostMeanSquare estimated posterior means of the squares of parameters
#' @slot gammaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot nuPostMeanSquare estimated posterior mean of the square of the parameter 
#' @slot omegaSquaredPostMeanSquare estimated posterior means of the squares of parameters
#' @slot rhoPostMeanSquare estimated posterior means of the squares of parameters
#' @slot sigmaSquaredPostMeanSquare estimated posterior means of the squares of parameters
#' @slot tauPostMeanSquare estimated posterior mean of the square of the parameter
#' @slot thetaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot xiPostMeanSquare posterior means of the squares of parameters
setClass("Chain",
  slots = list(
    psrf = "numeric",

    conjunctions = "integer",
    contrasts = "numeric",
    counts = "integer",
    design = "numeric",
    values = "numeric",

    diag = "character",
    ess = "integer",
    max_attempts = "numeric",
    nchains_diag = "integer",
    psrf_tol = "numeric",

    betas_update = "integer",
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

    C = "integer",
    countSums_g = "integer",
    countSums_n = "integer",
    designUnique = "numeric",
    designUniqueN = "integer",
    G = "integer",
    Greturn = "integer",
    GreturnEpsilon = "integer",
    J = "integer",
    L = "integer",
    Lupdate = "integer",
    N = "integer",
    Nreturn = "integer",
    NreturnEpsilon = "integer",
    probs = "numeric",
    seeds = "integer",

    a = "numeric",
    b = "numeric",
    c = "numeric",
    d = "numeric",
    k = "numeric",
    r = "numeric",
    s = "numeric",
    w = "numeric",

    beta = "numeric",
    epsilon = "numeric",
    gamma = "numeric",
    nu = "numeric",
    omegaSquared = "numeric",
    rho = "numeric",
    sigmaSquared = "numeric",
    tau = "numeric",
    theta = "numeric",
    xi = "numeric",

    betaStart = "numeric",
    epsilonStart = "numeric",
    gammaStart = "numeric",
    nuStart = "numeric",
    omegaSquaredStart = "numeric",
    rhoStart = "numeric",
    sigmaSquaredStart = "numeric",
    tauStart = "numeric",
    thetaStart = "numeric",
    xiStart = "numeric",

    betaPostMean = "numeric",
    epsilonPostMean = "numeric",
    gammaPostMean = "numeric",
    nuPostMean = "numeric",
    omegaSquaredPostMean = "numeric",
    rhoPostMean = "numeric",
    sigmaSquaredPostMean = "numeric",
    tauPostMean = "numeric",
    thetaPostMean = "numeric",
    xiPostMean = "numeric",

    betaPostMeanSquare = "numeric",
    epsilonPostMeanSquare = "numeric",
    gammaPostMeanSquare = "numeric",
    nuPostMeanSquare = "numeric",
    omegaSquaredPostMeanSquare = "numeric",
    rhoPostMeanSquare = "numeric",
    sigmaSquaredPostMeanSquare = "numeric",
    tauPostMeanSquare = "numeric",
    thetaPostMeanSquare = "numeric",
    xiPostMeanSquare = "numeric"
  )
)
