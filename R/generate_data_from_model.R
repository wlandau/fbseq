#' @include scenario.R
NULL

#' @title Function \code{generate_data_from_model}
#' @description Simulates count data fromt the Laplace-prior version of the model 
#' given hyperparameters. All hyperparameters must be given in the \code{truth} slot 
#' of the \code{Scenario} object argument.
#' @export
#' 
#' @return a list: a \code{Scenario} object with the scenario and a \code{Starts} object
#' with the true parameters that the count dataset in the \code{Scenario} object
#' was simulated from.
#' @param genes number of genes/features/response variables in the data
#' @param design Design matrix, with \code{N} rows and \code{L} columns.
#' Same format as in \code{edgeR::glmFit}, with RNA-seq libraries/samples as rows
#' and each \code{beta_{l, .}} as a column.
#' @param truth a \code{Starts} object containing the necessary initialization constants
#' and hyperparameters. All hyperparameters must be supplied.
#' @param priors character string, either a member of \code{special_beta_priors()} or 
#' another string for normal priors on the betas.
generate_data_from_model = function(genes, design, truth, priors = "Laplace"){
  libraries = nrow(design)

  stopifnot(libraries >= 3)
  stopifnot(libraries == nrow(design))

  truth@xi = numeric(0)
  truth@beta = numeric(0)
  for(l in 1:ncol(design)){
    
    if(priors == "Laplace"){
      xi = rexp(genes, rate = truth@k[1])
    } else if(priors == "t"){
      xi = 1/rgamma(genes, shape = truth@q[1], rate = truth@r[1])
    } else if(priors == "horseshoe"){
      xi = abs(rcauchy(genes))
    } else {
      xi = rep(1, genes)
    }
    xi[xi > 10] = sample(xi[xi <= 10], sum(xi > 10), replace = T)

    beta = rnorm(n = genes, mean = truth@theta[l], sd = sqrt(truth@sigmaSquared[l] * xi))
    truth@xi = c(truth@xi, xi)
    truth@beta = c(truth@beta, beta)
  }  

  truth@gamma = 1/rgamma(genes, shape = truth@nu/2, rate = truth@nu*truth@tau/2)

  truth@rho = rnorm(libraries, 0, sqrt(truth@omegaSquared))
  rhomat = matrix(rep(truth@rho, each = genes), ncol = libraries)

  gammat = matrix(rep(truth@gamma, times = libraries), ncol = libraries)
  epsilon = matrix(rnorm(libraries*genes, rhomat, sqrt(gammat)), ncol = libraries)

  beta = matrix(truth@beta, nrow = genes)
  lambda = t(design %*% t(beta))

  counts = matrix(rpois(genes*libraries, exp(lambda + epsilon)), nrow = genes)
  rownames(counts) = paste("gene_", 1:genes, sep="")
  colnames(counts) = rownames(design)

  truth@epsilon = as.vector(epsilon)
  Scenario(counts = counts, design = design, supplement = list(truth = truth))
}
