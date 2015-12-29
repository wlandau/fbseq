#' @include starts.R
NULL

#' @title Class \code{Chain}
#' @description the main storage object for a heterosis MCMC chain.
#' Create a new one with a call to \code{Chain()} and run MCMC
#' by feeding one into \code{single_chain()}.
#' @exportClass Chain
#'
#' @slot bound_names names of \code{bounds} slot in the original \code{Scenario} object
#' @slot contrast_names names of \code{contrasts} slot in the \code{Scenario} object
#' @slot gene_names names of genes, taken from the row names of the count data matrix
#' @slot library_names names of the libraries/samples, taken from the column names of the count data matrix
#' @slot proposition_names names of \code{propositions} slot in the \code{Scenario} object
#'
#' @slot bounds values to compare contrasts to. The comparison is to see if the contrast is greater than
#' its corresponding element in \code{bounds} (from the \code{Scenario} object)
#' @slot contrasts contrasts from the \code{Scenario} object.
#' @slot counts RNA-seq count data, flattened from a matrix
#' @slot design Design, flattened from the design matrix.
#' Original matrix must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' sets of gene-specific variables.
#' @slot propositions propositions of inequalities involving contrasts from the \code{Scenario} object.
#' @slot supplement a list containing supplementary information about the scenario: 
#' for example, how the data were simulated, if applicable
#'
#' @slot burnin MCMC burnin, the number of MCMC iterations to ignore at the beginning of each obj
#' @slot effects_update_beta bounds of l for which to update the beta_{l, g} parameters.
#' @slot genes_return Indices of genes whose parameter samples you want to return.
#' Applies to all gene-specific parameters except for the epsilons.
#' @slot genes_return_epsilon Indices of genes g for which epsilon_{n, g} is updated/returned.
#' @slot iterations Number of MCMC iterations after burnin for which selected parameter samples are kept.
#' Total MCMC iterations = burnin + thin * "iterations", and the whole "thin * iterations" portion
#' is used to calculate posterior means, mean squares, and probabilities.
#' @slot libraries_return Indices of RNA-seq libraries whose parameter samples you want to return.
#' Currently moot because there are no library-specific parameters other than the epsilons, but that
#' could change in future versions of the package.
#' @slot libraries_return_epsilon Indices of RNA-seq libraries n for which epsilon_{n, g} is updated/returned.
#' Applies to all library-specific parameters except for the epsilons.
#' @slot parameter_sets_return Character vector naming the variables whose MCMC samples 
#' you want to return
#' @slot parameter_sets_update Character vector naming the variables to calculate/update
#' during the MCMC.
#' @slot priors Names of the family of priors on the betas after integrating out the xi's. 
#' Can be any value returned by \code{special_beta_priors()}. All other bounds will default to the normal prior.
#' @slot thin MCMC thinning interval. \code{thin = 1} means parameter samples will be saved for every iterations
#' after burnin. \code{thin = 10} means parameter samples will be saved every 10th iteration after burnin.
#' Total MCMC iterations = burnin + thin * "iterations", and the whole "thin * iterations" portion
#' is used to calculate posterior means, mean squares, and probabilities.
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
#' @slot L number of columns in the original design matrix
#' @slot Lupdate_beta number of bounds of l for which to update the beta_{l, g} parameters.
#' @slot N number of libraries
#' @slot Nreturn number of libraries to return library-specific MCMC parameter samples for (except the epsilons)
#' @slot NreturnEpsilon number of libraries to return library-specific MCMC epsilon parameter samples
#' @slot P number of propositions involving contrasts
#' @slot probs estimated posterior probabilities of propositions in the \code{Scenario} object.
#' @slot seeds vector of N*G random number generator seeds
#' 
#' @slot a initialization constant  
#' @slot b initialization constant 
#' @slot c initialization constants 
#' @slot d initialization constant 
#' @slot h initialization constants 
#' @slot k initialization constants 
#' @slot q initialization constants 
#' @slot r initialization constants 
#' @slot s initialization constants 
#' 
#' @slot beta MCMC parameter samples
#' @slot epsilon MCMC parameter samples
#' @slot gamma MCMC parameter samples
#' @slot nu MCMC parameter samples
#' @slot sigmaSquared MCMC parameter samples
#' @slot tau MCMC parameter samples
#' @slot theta MCMC parameter samples
#' @slot xi MCMC parameter samples
#' 
#' @slot betaStart MCMC starting bounds
#' @slot epsilonStart MCMC starting bounds
#' @slot gammaStart MCMC starting bounds
#' @slot nuStart MCMC starting bounds
#' @slot sigmaSquaredStart MCMC starting bounds
#' @slot tauStart MCMC starting bounds
#' @slot thetaStart MCMC starting bounds
#' @slot xiStart MCMC starting bounds
#' 
#' @slot betaPostMean estimated posterior means
#' @slot epsilonPostMean estimated posterior means
#' @slot gammaPostMean estimated posterior means
#' @slot nuPostMean estimated posterior mean
#' @slot sigmaSquaredPostMean estimated posterior means
#' @slot tauPostMean estimated posterior mean
#' @slot thetaPostMean estimated posterior means
#' @slot xiPostMean estimated posterior means
#' 
#' @slot betaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot epsilonPostMeanSquare estimated posterior means of the squares of parameters
#' @slot gammaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot nuPostMeanSquare estimated posterior mean of the square of the parameter 
#' @slot sigmaSquaredPostMeanSquare estimated posterior means of the squares of parameters
#' @slot tauPostMeanSquare estimated posterior mean of the square of the parameter
#' @slot thetaPostMeanSquare estimated posterior means of the squares of parameters
#' @slot xiPostMeanSquare posterior means of the squares of parameters
#'
#' @slot betaSampler sampler option
#' @slot epsilonSampler sampler option
#' @slot gammaSampler sampler option
#' @slot nuSampler sampler option
#' @slot sigmaSquaredSampler sampler option
#' @slot tauSampler sampler option
#' @slot thetaSampler sampler option
#' @slot xiSampler sampler option
#' 
#' @slot betaTune tuning parameter
#' @slot epsilonTune tuning parameter
#' @slot gammaTune tuning parameter
#' @slot nuTune tuning parameter
#' @slot sigmaSquaredTune tuning parameter
#' @slot tauTune tuning parameter
#' @slot thetaTune tuning parameter
#' @slot xiTune tuning parameter
setClass("Chain",
  slots = list(
    bound_names = "character",
    contrast_names = "character",
    gene_names = "character",
    library_names = "character",
    proposition_names = "character",

    bounds = "numeric",
    contrasts = "numeric",
    counts = "integer",
    design = "numeric",
    propositions = "integer",
    supplement = "list",

    burnin = "integer",
    effects_update_beta = "integer",
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
    L = "integer",
    Lupdate_beta = "integer",
    N = "integer",
    Nreturn = "integer",
    NreturnEpsilon = "integer",
    P = "integer",
    probs = "numeric",
    seeds = "integer",

    a = "numeric",
    b = "numeric",
    c = "numeric",
    d = "numeric",
    h = "numeric",
    k = "numeric",
    q = "numeric",
    r = "numeric",
    s = "numeric",

    beta = "numeric",
    epsilon = "numeric",
    gamma = "numeric",
    nu = "numeric",
    sigmaSquared = "numeric",
    tau = "numeric",
    theta = "numeric",
    xi = "numeric",

    betaStart = "numeric",
    epsilonStart = "numeric",
    gammaStart = "numeric",
    nuStart = "numeric",
    sigmaSquaredStart = "numeric",
    tauStart = "numeric",
    thetaStart = "numeric",
    xiStart = "numeric",

    betaPostMean = "numeric",
    epsilonPostMean = "numeric",
    gammaPostMean = "numeric",
    nuPostMean = "numeric",
    sigmaSquaredPostMean = "numeric",
    tauPostMean = "numeric",
    thetaPostMean = "numeric",
    xiPostMean = "numeric",

    betaPostMeanSquare = "numeric",
    epsilonPostMeanSquare = "numeric",
    gammaPostMeanSquare = "numeric",
    nuPostMeanSquare = "numeric",
    sigmaSquaredPostMeanSquare = "numeric",
    tauPostMeanSquare = "numeric",
    thetaPostMeanSquare = "numeric",
    xiPostMeanSquare = "numeric",

    betaSampler = "integer",
    epsilonSampler = "integer",
    gammaSampler = "integer",
    nuSampler = "integer",
    sigmaSquaredSampler = "integer",
    tauSampler = "integer",
    thetaSampler = "integer",
    xiSampler = "integer",

    betaTune = "numeric",
    epsilonTune = "numeric",
    gammaTune = "numeric",
    nuTune = "numeric",
    sigmaSquaredTune = "numeric",
    tauTune = "numeric",
    thetaTune = "numeric",
    xiTune = "numeric"
  )
)
