# fbseq

`fbseq` is an R package for analyzing bowtie-preprocessed RNA-sequencing datasets with arbitrary experimental designs and computing arbitrary posterior probabilities of interest. The computation is accelerated with CUDA. For an overview of the intended applications, the model, and the methodology, see the [model vignette](https://github.com/wlandau/fbseq/blob/master/vignettes/model.html). For installation and usage instructions, see the [tutorial vignette](https://github.com/wlandau/fbseq/blob/master/vignettes/tutorial.html).For best viewing of either html vignette, download it to your desktop and open it with a browser.

# Read the [model vignette](https://github.com/wlandau/fbseq/blob/master/vignettes/model.html) first. 

An understanding of the underlying hierarchical model is important for understanding how to use this package. For best viewing, download the html file to your desktop and then open it with a browser.

# Check your system.

You need to have at least R $\ge$ 3.2.0, along with the R packages `coda`, `ggplot2`, `methods`, `reshape2`, and `knitr`. All are available through the [Comprehensive R Archive Network (CRAN](https://cran.r-project.org/). With those requirements met, you can install `fbseq`, load it in an R session, create input, and analyze output. To actually run the underlying Markov chain Monte Carlo (MCMC) procedure, however, you need access to a machine with a [CUDA-enabled GPU](http://www.nvidia.com/object/cuda_home_new.html), along with the [`fbseqCUDA` package](https://github.com/wlandau/fbseqCUDA). `fbseqCUDA` is the internal engine of `fbseq`, and it is implemented in CUDA to provide necessary acceleration for the MCMC. `fbseq` and `fbseqCUDA` are kept separate for convenience. For example, you can set up input with `fbseq` and without `fbseqCUDA` on a low-end laptop, run the main algorithm remotely with both `fbseq` and `fbseqCUDA` on a CUDA-enabled [cloud computing enterprise](http://www.nvidia.com/object/gpu-cloud-computing-services.html) such as [Amazon Web Services](http://aws.amazon.com/ec2/instance-types/), and analyze the output locally with `fbseq` and without `fbseqCUDA`. See the [`fbseqCUDA` installation vignette](https://github.com/wlandau/fbseqCUDA/blob/master/vignettes/install.html) for more details. For best viewing, download the html vignette to your desktop and then open it with a browser.

# Install `fbseq`

## Option 1: install a stable release (recommended).

For `fbseq`, you can navigate to a [list of stable releases](https://github.com/wlandau/fbseq/releases) on the project's [GitHub page](https://github.com/wlandau/fbseq). Download the desired `tar.gz` bundle, then install it either with `install.packages(..., repos = NULL, type="source")` from within R  `R CMD INSTALL` from the Unix/Linux command line.

## Option 2: use `install_github` to install the development version.

For this option, you need the `devtools` package, available from CRAN or GitHub. Open R and run 

```{r, eval=F}
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

# Install `fbseqCUDA` 

You can only install `fbseqCUDA` on a machine with CUDA and a CUDA-capable graphics processing unit (GPU). The `fbseqCUDA` package is necessary for running the MCMC, but not for preparing input or analyzing output. For installation instructions, see the [`fbseqCUDA` installation vignette](https://github.com/wlandau/fbseqCUDA/blob/master/vignettes/install.html).
