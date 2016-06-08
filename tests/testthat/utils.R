data(paschold)
paschold@counts = paschold@counts[1:20,]
cf = Configs(iterations = 50, burnin = 50, thin = 1, genes_return = 1:20, verbose = 0)
ch = Chain(paschold, cf)
contrasts = paschold@contrasts
propositions = paschold@propositions
threads = 1

skip_if_missing_backend = function(backend){
  if(!requireNamespace(backends()[backend], quietly = T))
    skip(paste("Backend", backend, "missing."))
}
