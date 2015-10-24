#' @title Class \code{Chain}
#' @description the main storage object for a heterosis MCMC chain.
#' Create a new one with a call to \code{Chain()} and run MCMC
#' by feeding one into \code{single_chain()}.
#' @exportClass Chain
#'
#' @slot diag convergence diagnostic to use. Can be "gelman" or "none".
#' @slot ess Minimum effective sample size for all parameters
#' @slot psrf_tol upper threshold for Gelman-Rubin potential scale reduction factors (if diag is "gelman")
#' @slot max_attempts Maximum number of retries for assessing convergence and generating enough effective samples.
#' Can be set to Inf to run indefinitely.
#' @slot nchains_diag number of independent chains to run (including this one) to use
#' convergence diagnostics that require multiple chains. 
#' @slot psrf Gelman-Rubin potential scale reduction factors, if running \code{run_gelman_mcmc()}
#' 
#' @slot counts vector of RNA-seq read counts. Type \code{matrix(chain@@counts, ncol = chain@@N)}
#' to recover the original RNA-seq count data matrix with features/genes
#' as rows and samples/libraries as columns with .
#' @slot group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for parent 2,
#' and 3 for the hybrid.
#' 
#' @slot sums_n Sample/library-specific sums of counts.
#' @slot sums_g Feature/gene-specific sums of counts.
#' 
#' @slot returns named integer vector of 0/1 values naming
#' the variables whose MCMC samples you want to return.
#' @slot updates named integer vector of 0/1 values naming 
#' the variables to compute/update during the MCMC.
#' 
#' @slot samples_return indices of RNA-seq samples/libraries whose parameter samples you want to return.
#' @slot features_return indices of features/genes whose parameter samples you want to return.
#'
#' @slot samples_return_eps Indices of RNA-seq samples/libraries n for which epsilon_{n, g} is updated/returned.
#' @slot features_return_eps Indices of features/genes g for which epsilon_{n, g} is updated/returned.
#' 
#' @slot M number of MCMC iterations (not including burnin or thinning)
#' @slot N number of RNA-seq samples/libraries
#' @slot G number of features/genes
#'
#' @slot Nreturn length of slot \code{samples_return}.
#' @slot Greturn length of slot \code{features_return}.
#'
#' @slot NreturnEps length of slot \code{samples_return_eps}.
#' @slot GreturnEps length of slot \code{features_return_eps}.
#'
#' @slot burnin MCMC burnin, the number of MCMC iterations to ignore at the beginning of each chain
#' @slot mphtol tolerance theshold for mid-parent heterosis (on a log base e scale)
#' @slot thin MCMC thinning interval, number of iterations to skip in between iterations to return.
#' @slot verbose number of times to print out progress during burnin and the actual MCMC.
#' If \code{verbose} > 0, then progress messages will also print during setup and cleanup.
#' @slot seeds Random number generator seeds.
#'
#' @slot hph gene/feature-specific high-parent heterosis probabilities. The order of the 
#' probabilities respects the ordering of features/genes in the original RNA-seq dataset.
#' @slot lph gene/feature-specific low-parent heterosis probabilities. The order of the 
#' probabilities respects the ordering of features/genes in the original RNA-seq dataset.
#' @slot mph gene/feature-specific mid-parent heterosis probabilities. The order of the 
#' probabilities respects the ordering of features/genes in the original RNA-seq dataset.
#' 
#' @slot phiPrior Integer indicating the family of priors on the phi_g's after integrating out the xi_phi's. 
#' The integer is either 0 for normal or the index of the corresponding prior in the 
#' return value of \code{alternate_priors()}.
#' @slot alpPrior Integer indicating the family of priors on the alp_g's after integrating out the xi_alp's. 
#' The integer is either 0 for normal or the index of the corresponding prior in the 
#' return value of \code{alternate_priors()}.
#' @slot delPrior Integer indicating the family of priors on the del_g's after integrating out the xi_del's. 
#' The integer is either 0 for normal or the index of the corresponding prior in the 
#' return value of \code{alternate_priors()}.
#'
#' @slot dRho initialization constant
#' @slot aRho initialization constant
#' @slot bRho initialization constant
#' @slot dGam initialization constant
#' @slot aGam initialization constant
#' @slot bGam initialization constant
#' @slot cPhi initialization constant
#' @slot cAlp initialization constant
#' @slot cDel initialization constant
#' @slot sPhi initialization constant
#' @slot sAlp initialization constant
#' @slot sDel initialization constant
#' @slot kPhi initialization constant
#' @slot kAlp initialization constant
#' @slot kDel initialization constant
#' @slot rPhi initialization constant
#' @slot rAlp initialization constant
#' @slot rDel initialization constant
#' 
#' @slot nuRho vector of model parameter samples
#' @slot nuGam vector of model parameter samples
#' @slot tauRho vector of model parameter samples
#' @slot tauGam vector of model parameter samples
#' @slot thePhi vector of model parameter samples
#' @slot theAlp vector of model parameter samples
#' @slot theDel vector of model parameter samples
#' @slot sigPhi vector of model parameter samples
#' @slot sigAlp vector of model parameter samples
#' @slot sigDel vector of model parameter samples
#' 
#' @slot phi vector of model parameter samples
#' @slot alp vector of model parameter samples
#' @slot del vector of model parameter samples
#' @slot rho vector of model parameter samples
#' @slot gam vector of model parameter samples
#' @slot xiPhi vector of model parameter samples
#' @slot xiAlp vector of model parameter samples
#' @slot xiDel vector of model parameter samples
#' @slot eps vector of model parameter samples (flattened from a 2D array)
#' 
#' @slot nuRhoStart vector of starting values
#' @slot nuGamStart vector of starting values
#' @slot tauRhoStart vector of starting values
#' @slot tauGamStart vector of starting values
#' @slot thePhiStart vector of starting values
#' @slot theAlpStart vector of starting values
#' @slot theDelStart vector of starting values
#' @slot sigPhiStart vector of starting values
#' @slot sigAlpStart vector of starting values
#' @slot sigDelStart vector of starting values
#'
#' @slot phiStart vector of starting values
#' @slot alpStart vector of starting values
#' @slot delStart vector of starting values
#' @slot rhoStart vector of starting values
#' @slot gamStart vector of starting values
#' @slot xiPhiStart vector of starting values
#' @slot xiAlpStart vector of starting values
#' @slot xiDelStart vector of starting values
#' @slot epsStart vector of starting values (flattened from a 2D array)
#'
#' @slot nuRhoPostMean vector of posterior means after burnin
#' @slot nuGamPostMean vector of posterior means after burnin
#' @slot tauRhoPostMean vector of posterior means after burnin
#' @slot tauGamPostMean vector of posterior means after burnin
#' @slot thePhiPostMean vector of posterior means after burnin
#' @slot theAlpPostMean vector of posterior means after burnin
#' @slot theDelPostMean vector of posterior means after burnin
#' @slot sigPhiPostMean vector of posterior means after burnin
#' @slot sigAlpPostMean vector of posterior means after burnin
#' @slot sigDelPostMean vector of posterior means after burnin
#'
#' @slot phiPostMean vector of posterior means after burnin
#' @slot alpPostMean vector of posterior means after burnin
#' @slot delPostMean vector of posterior means after burnin
#' @slot rhoPostMean vector of posterior means after burnin
#' @slot gamPostMean vector of posterior means after burnin
#' @slot xiPhiPostMean vector of posterior means after burnin
#' @slot xiAlpPostMean vector of posterior means after burnin
#' @slot xiDelPostMean vector of posterior means after burnin
#' @slot epsPostMean vector of posterior means after burnin (flattened from a 2D array)
#'
#' @slot nuRhoPostMeanSq vector of posterior means of squares after burnin
#' @slot nuGamPostMeanSq vector of posterior means of squares after burnin
#' @slot tauRhoPostMeanSq vector of posterior means of squares after burnin
#' @slot tauGamPostMeanSq vector of posterior means of squares after burnin
#' @slot thePhiPostMeanSq vector of posterior means of squares after burnin
#' @slot theAlpPostMeanSq vector of posterior means of squares after burnin
#' @slot theDelPostMeanSq vector of posterior means of squares after burnin
#' @slot sigPhiPostMeanSq vector of posterior means of squares after burnin
#' @slot sigAlpPostMeanSq vector of posterior means of squares after burnin
#' @slot sigDelPostMeanSq vector of posterior means of squares after burnin
#'
#' @slot phiPostMeanSq vector of posterior means of squares after burnin
#' @slot alpPostMeanSq vector of posterior means of squares after burnin
#' @slot delPostMeanSq vector of posterior means of squares after burnin
#' @slot rhoPostMeanSq vector of posterior means of squares after burnin
#' @slot gamPostMeanSq vector of posterior means of squares after burnin
#' @slot xiPhiPostMeanSq vector of posterior means of squares after burnin
#' @slot xiAlpPostMeanSq vector of posterior means of squares after burnin
#' @slot xiDelPostMeanSq vector of posterior means of squares after burnin
#' @slot epsPostMeanSq vector of posterior means of squares after burnin (flattened from a 2D array)
setClass("Chain",
  slots = list(
    diag = "character",
    ess = "integer",
    psrf_tol = "numeric",
    max_attempts = "numeric",
    nchains_diag = "integer",
    psrf = "numeric",

    counts = "integer",
    group = "integer",

    sums_n = "integer",
    sums_g = "integer",

    returns = "integer",
    updates = "integer",

    samples_return = "integer",
    features_return = "integer",

    samples_return_eps = "integer",
    features_return_eps = "integer",

    M = "integer",
    N = "integer",
    G = "integer",

    Nreturn = "integer",
    Greturn = "integer",

    NreturnEps = "integer",
    GreturnEps = "integer",

    burnin = "integer",
    thin = "integer",
    verbose = "integer",
    seeds = "integer",
    mphtol = "numeric",

    hph = "numeric",
    lph = "numeric",
    mph = "numeric",

    phiPrior = "integer",
    alpPrior = "integer",
    delPrior = "integer",

    dRho = "numeric",
    dGam = "numeric",
    aRho = "numeric",
    aGam = "numeric",
    bRho = "numeric",
    bGam = "numeric",
    cPhi = "numeric",
    cAlp = "numeric",
    cDel = "numeric",
    sPhi = "numeric",
    sAlp = "numeric",
    sDel = "numeric",
    kPhi = "numeric",
    kAlp = "numeric",
    kDel = "numeric",
    rPhi = "numeric",
    rAlp = "numeric",
    rDel = "numeric",

    nuRho = "numeric",
    nuGam = "numeric",
    tauRho = "numeric",
    tauGam = "numeric",
    thePhi = "numeric",
    theAlp = "numeric",
    theDel = "numeric",
    sigPhi = "numeric",
    sigAlp = "numeric",
    sigDel = "numeric",

    phi = "numeric",
    alp = "numeric",
    del = "numeric",  
    rho = "numeric",
    gam = "numeric",
    xiPhi = "numeric",
    xiAlp = "numeric",
    xiDel = "numeric",
    eps = "numeric",

    nuRhoStart = "numeric",
    nuGamStart = "numeric",
    tauRhoStart = "numeric",
    tauGamStart = "numeric",
    thePhiStart = "numeric",
    theAlpStart = "numeric",
    theDelStart = "numeric",
    sigPhiStart = "numeric",
    sigAlpStart = "numeric",
    sigDelStart = "numeric",

    phiStart = "numeric",
    alpStart = "numeric",
    delStart = "numeric",
    rhoStart = "numeric",
    gamStart = "numeric",
    xiPhiStart = "numeric",
    xiAlpStart = "numeric",
    xiDelStart = "numeric",
    epsStart = "numeric",

    nuRhoPostMean = "numeric",
    nuGamPostMean = "numeric",
    tauRhoPostMean = "numeric",
    tauGamPostMean = "numeric",
    thePhiPostMean = "numeric",
    theAlpPostMean = "numeric",
    theDelPostMean = "numeric",
    sigPhiPostMean = "numeric",
    sigAlpPostMean = "numeric",
    sigDelPostMean = "numeric",

    phiPostMean = "numeric",
    alpPostMean = "numeric",
    delPostMean = "numeric",
    rhoPostMean = "numeric",
    gamPostMean = "numeric",
    xiPhiPostMean = "numeric",
    xiAlpPostMean = "numeric",
    xiDelPostMean = "numeric",
    epsPostMean = "numeric", 

    nuRhoPostMeanSq = "numeric",
    nuGamPostMeanSq = "numeric",
    tauRhoPostMeanSq = "numeric",
    tauGamPostMeanSq = "numeric",
    thePhiPostMeanSq = "numeric",
    theAlpPostMeanSq = "numeric",
    theDelPostMeanSq = "numeric",
    sigPhiPostMeanSq = "numeric",
    sigAlpPostMeanSq = "numeric",
    sigDelPostMeanSq = "numeric",

    phiPostMeanSq = "numeric",
    alpPostMeanSq = "numeric",
    delPostMeanSq = "numeric",
    rhoPostMeanSq = "numeric",
    gamPostMeanSq = "numeric",
    xiPhiPostMeanSq = "numeric",
    xiAlpPostMeanSq = "numeric",
    xiDelPostMeanSq = "numeric",
    epsPostMeanSq = "numeric"
  )
)
