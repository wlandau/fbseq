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
#' @param design Gene-specific design matrix. Must contain only 0's, 1's, and -1's.
#' Must have rows corresponding to colums/libraries in RNA-seq data and colums corresponding to
#' gene-specific variables.
#' Can be among "Laplace", "t", or "horseshoe". All other values will default to the normal prior.
generate_data = function(libraries = 12, genes = 3.5e4, 
          design = cbind(rep(1, libraries), rep(c(1, -1, 1), each = floor(libraries/3)), rep(c(-1, 1, 1), each = floor(libraries/3))),
          starts = Starts(nuGamma = 5, nuRho = 5, omega = c(1, 0.5, 0.5), tauGamma = 1, tauRho = 0.1, theta = c(3, 0, 0))){

  stopifnot(libraries >= 3)

  starts@xi = rep(1, ncol(design)*genes)
  for(l in 1:ncol(design))
    starts@beta = c(starts@beta, rnorm(genes, starts@theta[l], sqrt(starts@omega[l])))
  
  starts@gamma = 1/rgamma(genes, shape = starts@nuGam/2, 
    rate = starts@nuGam*starts@tauGam^2/2)
  starts@rho = 1/rgamma(libraries, shape = starts@nuRho/2, 
    rate = starts@nuRho*starts@tauRho^2/2)

  rhomat = matrix(rep(starts@rho, each = genes), ncol = libraries)
  gammat = matrix(rep(starts@gamma, times = libraries), ncol = libraries)
  epsilon = matrix(rnorm(libraries*genes, 0, sqrt(rhomat*gammat)), ncol = libraries)

  beta = matrix(starts@beta, nrow = genes)
  lambda = t(design %*% t(beta))

  counts = matrix(rpois(genes*libraries, exp(lambda + epsilon)), nrow = genes)
  rownames(counts) = paste("gene_", 1:genes, sep="")
  colnames(counts) = paste("library_", 1:libraries, sep="")

  starts@epsilon = as.vector(epsilon)
  return(list(counts = counts, group = group, truth = starts))
}
