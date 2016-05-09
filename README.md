# Introduction

The [`fbseq`](https://github.com/wlandau/fbseq)  package part of a collection of packages for the fully Bayesian analysis of RNA-sequencing count data, where a hierarchical model is fit with Markov chain Monte Carlo (MCMC). [`fbseq`](https://github.com/wlandau/fbseq)  is the user interface, and it contains top level functions for calling the MCMC and analyzing output. The other packages, [`fbseqSerial`](https://github.com/wlandau/fbseqSerial)  and [`fbseqCUDA`](https://github.com/wlandau/fbseqCUDA) , are backend packages that run the MCMC behind the scenes. Only one is required.  [`fbseqSerial`](https://github.com/wlandau/fbseqSerial)  can run on most machines, but it is slow for large datasets. [`fbseqCUDA`](https://github.com/wlandau/fbseqCUDA)  requires special hardware (i.e. a CUDA general-purpose graphics processing unit), but it's much faster due to parallel computing. 

# Read the [model vignette](https://github.com/wlandau/fbseq/blob/master/vignettes/model.html) first. 

An understanding of the underlying hierarchical model is important for understanding how to use this package. For best viewing, download the html file to your desktop and then open it with a browser.

# Check your system.

See the "Depends", "Imports", and "Suggests" fields of the "package's [DESCRIPTION](https://github.com/wlandau/fbseq/blob/master/DESCRIPTION) R version and R package requirements.

# Install [`fbseq`](https://github.com/wlandau/fbseq) 

## Option 1: install a stable release (recommended).

For [`fbseq`](https://github.com/wlandau/fbseq) , you can navigate to a [list of stable releases](https://github.com/wlandau/fbseq/releases) on the project's [GitHub page](https://github.com/wlandau/fbseq). Download the desired `tar.gz` bundle, then install it either with `install.packages(..., repos = NULL, type="source")` from within R  `R CMD INSTALL` from the Unix/Linux command line.

## Option 2: use `install_github` to install the development version.

For this option, you need the `devtools` package, available from CRAN or GitHub. Open R and run 

```
library(devtools)
install_github("wlandau/fbseq")
```

## Option 3: build the development version from the source.

Open a command line program such as Terminal in Mac/Linux and enter the following commands.

```
git clone git@github.com:wlandau/fbseq.git
R CMD build fbseq
R CMD INSTALL ...
```

where `...` is replaced by the name of the tarball produced by `R CMD build`. 

# Install an MCMC backend package

[`fbseqSerial`](https://github.com/wlandau/fbseqSerial)  and [`fbseqCUDA`](https://github.com/wlandau/fbseqCUDA) , are backend packages that run the MCMC behind the scenes. Only one is required.  [`fbseqSerial`](https://github.com/wlandau/fbseqSerial)  can run on most machines, but it is slow for large datasets. [`fbseqCUDA`](https://github.com/wlandau/fbseqCUDA)  requires special hardware (i.e. a CUDA general-purpose graphics processing unit), but it's much faster due to parallel computing. Installation is similar to that of [`fbseq`](https://github.com/wlandau/fbseq)  and is detailed in their respective `README.md` files and package vignettes.

# Quick start

After installing [`fbseq`](https://github.com/wlandau/fbseq)  and [`fbseqSerial`](https://github.com/wlandau/fbseqSerial) , the following should take a couple seconds to run. The example walks through an example data analysis and shows a few key features of the package. For more specific operational details, see the [tutorial vignette](https://github.com/wlandau/fbseq/blob/master/vignettes/tutorial.html).

```
library(fbseq)

back_end = "serial" # change this to "CUDA" to use fbseqCUDA as the backend

# Example RNA-seq dataset wrapped in an S4 class.
data(paschold) 

# Use a small subset of the data (few genes).
paschold@counts = head(paschold@counts, 20) 

# Run the MCMC for only a few iterations.
configs = Configs(burnin = 10, iterations = 10, thin = 1) 

# Set up the MCMC.
chain = Chain(paschold, configs) 

# Run 4 dispersed independent Markov chains.
chain_list = fbseq(chain, backend = back_end)

# Monitor convergence with in all parameters with 
# Gelman-Rubin potential scale reduction factors
head(psrf(chain_list)) 

# Means, standard deviations, and credible intervals of all parameters 
head(estimates(chain_list))

# Gene-specific (row-specific) posterior probabilities
head(probs(chain_list))

# Monte Carlo samples on a small subset of parameters
m = mcmc_samples(chain_list) 
dim(m)
m[1:5, 1:5]

# Set max_iter=Inf below to automatically run until convergence.
max_iter = 3
iter = 1
chain@verbose = as.integer(0) # turn off console messages
chain_list = fbseq(chain, backend = back_end)
while(iter < max_iter){
  chain_list = lapply(chain_list, fbseq, additional_chains = 0, backend = back_end)
  gelman = psrf(chain_list)
  # saveRDS(chain_list, "chain_list_so_far.rds") # option to save progress
  if(all(gelman < 1.1)) break
  iter = iter + 1
}
```
