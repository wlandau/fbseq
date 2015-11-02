#' @title Function \code{generate_data}
#' @description Simulates heterosis data from prespecified hyperparameters.
#' @export
#' 
#' @return a list
#' \describe{
#'   \item{\code{counts}}{The generated RNA-seq counts matrix.}
#'   \item{\code{design}}{The design matrix for gene-specific effects.}
#'   \item{\code{truth}}{A \code{Starts} object containing the values used to simulate the data}
#' }
#' @param starts A \code{Starts} object containing the necessary initialization constants
#' and hyperparameters.
#' @param libraries number of libraries/libraries in the data
#' @param genes number of genes/genes in the data
#' @param design Gene-specific design matrix.
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
#' Can be among "Laplace", "t", or "horseshoe". All other values will default to the normal prior.
generate_data = function(libraries = 12, genes = 3.5e4, 
          design = cbind(rep(1, libraries), rep(c(1, -1, 1), each = floor(libraries/3)), rep(c(-1, 1, 1), each = floor(libraries/3))),
          starts = Starts(nu = 10, omegaSquared = 0.01, pi = c(0.5, 1, 1), sigmaSquared = c(2.25, 0.25, 0.25), tau = 0.1, 
          theta = c(5, 0, 0))){

  stopifnot(libraries >= 3)
  starts@xi = rep(1, ncol(design)*genes)
  starts@p = c(1, 0, 0)

  starts@beta = starts@delta = numeric(0)
  for(l in 1:ncol(design)){
    deltas = sample(0:1, size = genes, replace = T, prob = c(1 - starts@pi[l], starts@pi[l]))
    starts@beta = c(starts@beta, rnorm(genes, deltas * starts@theta[l], sqrt(starts@sigmaSquared[l])))
    starts@delta = c(starts@delta, deltas)
  }

  starts@gamma = 1/rgamma(genes, shape = starts@nu/2, 
    rate = starts@nu*starts@tau/2)

  starts@rho = rnorm(libraries, 0, sqrt(starts@omegaSquared))
  rhomat = matrix(rep(starts@rho, each = genes), ncol = libraries)

  gammat = matrix(rep(starts@gamma, times = libraries), ncol = libraries)
  epsilon = matrix(rnorm(libraries*genes, rhomat, sqrt(gammat)), ncol = libraries)

  beta = matrix(starts@beta, nrow = genes)
  lambda = t(design %*% t(beta))

  counts = matrix(rpois(genes*libraries, exp(lambda + epsilon)), nrow = genes)
  rownames(counts) = paste("gene_", 1:genes, sep="")
  colnames(counts) = paste("library_", 1:libraries, sep="")

  starts@epsilon = as.vector(epsilon)
  return(list(counts = counts, design = design, truth = starts))
}
