# library(fbseq); library(testthat); backend = backend
source("utils.R")

for(backend in c("OpenMP", "CUDA")){

context(paste(backend, "horseshoe"))
threads = 1 + (fbseqOpenMP::OpenMP_working() & (backend == "OpenMP"))

test_that(paste(backend, "thetas are sampled iff priors != \"horseshoe\"."), {
  skip_if_missing_backend(backend)
  cf@priors = c("horseshoe", "Laplace", "horseshoe", "t", "normal")
  cf@parameter_sets_update = cf@parameter_sets_return = parameters()
  cf@theta_update = which(cf@priors != "horseshoe")
  ch = Chain(paschold, cf)
  ch@thetaStart[2:ncol(paschold@design)] = 0
  obj = fbseq(ch, backend = backend, additional_chains = 3, threads = threads)
  m = mcmc_samples(obj)
  theta = m[,grep("theta", colnames(m))]
  p = psrf(obj)
  v = unname(apply(theta, 2, var))
  expect_equal(v[cf@priors == "horseshoe"], rep(0, sum(cf@priors == "horseshoe")))
  expect_false(isTRUE(all.equal(v[cf@priors != "horseshoe"], 
    rep(0, sum(cf@priors != "horseshoe")))))
  for(i in 1:length(cf@priors)){
    xi = m[,grep(paste0("xi_", i), colnames(m))]
    u = unique(as.numeric(as.matrix(xi)))
    if(cf@priors[i] == "normal"){
      expect_equal(u, 1)
    } else {
      expect_true(length(u) > 1)
      expect_false(all(u == 1))
    }
  }
})
}
