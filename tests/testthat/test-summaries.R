# library(fbseq); library(testthat); backend = "OpenMP"
source("utils.R")

for(backend in c("OpenMP", "CUDA")){

context(paste(backend, "summaries"))
if(backend == "OpenMP") threads = 1 + fbseqOpenMP::OpenMP_working()

test_that(paste(backend, "contrast means are calculated correctly."), {
  skip_if_missing_backend(backend)
  obj = fbseq(ch, backend = "OpenMP", additional_chains = 3, threads = threads)
  samples = Reduce("+", lapply(obj, mcmc_samples))/4
  contrast_est = contrast_estimates(obj)
  for(g in 1:dim(paschold@counts)[1]){ 
    beta = as.matrix(samples[,grep(paste0("beta_[1-5]_", g, "$"), colnames(samples))])
    for(c in 1:length(contrasts)){
      contrast_g = beta %*% contrasts[[c]]
      n = paste0(names(contrasts)[c], "_", rownames(paschold@counts)[g])
      expect_equal(contrast_est[n, "mean"], mean(contrast_g))
    }
  }
})

test_that(paste(backend, "contrast mean squares are calculated correctly."), {
  skip_if_missing_backend(backend)
  obj = fbseq(ch, backend = "OpenMP", additional_chains = 0, threads = threads)
  samples = mcmc_samples(obj)
  contrast_est = contrast_estimates(obj)
  for(g in 1:dim(paschold@counts)[1]){ 
    beta = as.matrix(samples[,grep(paste0("beta_[1-5]_", g, "$"), colnames(samples))])
    for(c in 1:length(contrasts)){
      contrast_g = beta %*% contrasts[[c]]
      n = paste0(names(contrasts)[c], "_", rownames(paschold@counts)[g])
      expect_equal(contrast_est[n, "mean"], mean(contrast_g))
      expect_equal(contrast_est[n, "sd"], sd(contrast_g))
    }
  }
})

test_that(paste(backend, "functions estimates and contrast_estimates works with multiple chains."), {
  skip_if_missing_backend(backend)
  obj = fbseq(ch, backend = "OpenMP", additional_chains = 3, threads = threads)
  for(f in c("estimates", "contrast_estimates")){
    est = get(f)(obj)
    ests = lapply(obj, get(f))
    means = lapply(ests, function(x) x$mean)
    expect_equal(est$mean, Reduce("+", means)/4)
  }
})

test_that(paste(backend, "parameter sample means and mean squares are calculated correctly."), {
  skip_if_missing_backend(backend)
  cf = Configs(iterations = 50, burnin = 50, thin = 1, 
    genes_return = 1:20, libraries_return = 1:16, 
    genes_return_epsilon = 1:20, libraries_return_epsilon = 1:16,
    verbose = 0, priors = "Laplace")
  ch = Chain(paschold, cf)
  obj = fbseq(ch, backend = "OpenMP", additional_chains = 0, threads = threads)
  est = estimates(obj)
  samples = mcmc_samples(obj)
  n = colnames(samples)
  expect_equal(est[n,"mean"], unname(apply(samples, 2, mean)[n]))
  expect_equal(est[n,"sd"], unname(apply(samples, 2, sd)[n]))
})

}
