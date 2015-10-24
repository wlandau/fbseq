#' @title Function \code{generate_data}
#' @description Simulates heterosis data from prespecified hyperparameters.
#' Intended for testing and demonstration only.
#' Default hyperparameters are available as default arguments.
#' @details Each feature/gene must have at least one nonzero count per genetic variety.
#' The offending zero-count genes are replaced with bootstrap samples of the other genes.
#' @export
#' 
#' @return a list
#' \describe{
#'   \item{\code{counts}}{The generated RNA-seq counts matrix.}
#'   \item{\code{group}}{The group structure.}
#'   \item{\code{truth}}{A \code{Starts} object containing the values used to simulate the data}
#' }
#' @param starts A \code{Starts} object containing the necessary initialization constants
#' and hyperparameters.
#' @param samples number of samples/libraries in the data
#' @param features number of features/genes in the data
#' @param group Experimental design. A vector of integers,
#' one for each RNA-seq sample/library, denoting the genetic
#' variety of that sample. You must use 1 for parent 1, 2 for the hybrid,
#' and 3 for parent 2.
#' @param phiPrior Name of the family of priors on the phi_g's after integrating out the xi_phi's. 
#' Can be "Laplace", "t", or "horseshoe". All other values will default to the normal prior.
#' @param alpPrior Name of the family of priors on the alp_g's after integrating out the xi_alp's. 
#' Can be "Laplace", "t", or "horseshoe". All other values will default to the normal prior.
#' @param delPrior Name of the family of priors on the del_g's after integrating out the xi_del's. 
#' Can be "Laplace", "t", or "horseshoe". All other values will default to the normal prior.
generate_data = function(samples = 12, features = 3.5e4, 
          group = c(rep(1:3, each = floor(samples/3)), rep(3, samples %% 3)),
          starts = Starts(nuRho = 5, nuGam = 5, tauRho = 0.1, thePhi = 3, theAlp = 0, theDel = 0,
          sigPhi = 1, sigAlp = 0.5, sigDel = 0.5), phiPrior = "normal", alpPrior = "normal", delPrior = "normal"){

  stopifnot(samples >= 3)

  for(p in c("phi", "alp", "del")){
    prior = get(paste0(p, "Prior"))
    cap = paste0(toupper(substring(p, 1, 1)), substring(p, 2),collapse="")
    k = slot(starts, paste0("k", cap))
    r = slot(starts, paste0("r", cap))

    if(prior == "Laplace"){
      x = rexp(features, rate = k)
    } else if (prior == "t") {
      x = 1/rgamma(features, shape = k, rate = r)
    } else if (prior == "horseshoe") {
      x = abs(rcauchy(2*features, location = 0, scale = 1))
      x = sample(x[x < 5], features, replace = T)
    } else {
      x = rep(1, features)
    }

    slot(starts, paste0("xi", cap)) = x
  }

  starts@phi = rnorm(features, starts@thePhi, starts@sigPhi * sqrt(starts@xiPhi)) 
  starts@alp = rnorm(features, starts@theAlp, starts@sigAlp * sqrt(starts@xiAlp)) 
  starts@del = rnorm(features, starts@theDel, starts@sigDel * sqrt(starts@xiDel))
  starts@rho = sqrt(1/rgamma(samples, shape = starts@nuRho/2, 
    rate = starts@nuRho*starts@tauRho^2/2))
  starts@gam = sqrt(1/rgamma(features, shape = starts@nuGam/2, 
    rate = starts@nuGam*starts@tauGam^2/2))

  rhomat = matrix(rep(starts@rho, each = features), ncol = samples)
  gammat = matrix(rep(starts@gam, times = samples), ncol = samples)
  eps = matrix(rnorm(samples*features, 0, rhomat*gammat), ncol = samples)

  lambda = matrix(0, nrow = features, ncol = samples)
  lambda = sweep(lambda, 1, starts@phi, "+")

  n = table(group)[as.character(1:3)]

  lambda[, group == 1] = sweep(matrix(lambda[, group == 1], ncol = n[1]), 1, starts@alp, "+")
  lambda[, group == 2] = sweep(matrix(lambda[, group == 2], ncol = n[2]), 1, starts@alp, "-")
  lambda[, group == 3] = sweep(matrix(lambda[, group == 3], ncol = n[3]), 1, starts@alp, "+")

  lambda[, group == 1] = sweep(matrix(lambda[, group == 1], ncol = n[1]), 1, starts@del, "-")
  lambda[, group == 2] = sweep(matrix(lambda[, group == 2], ncol = n[2]), 1, starts@del, "+")
  lambda[, group == 3] = sweep(matrix(lambda[, group == 3], ncol = n[3]), 1, starts@del, "+")

  counts = matrix(rpois(features*samples, exp(lambda + eps)), nrow = features)
  rownames(counts) = paste("feature_", 1:features, sep="")
  colnames(counts) = paste("sample_", 1:samples, sep="")

  starts@eps = as.vector(eps)
  return(list(counts = counts, group = group, truth = starts))
}
