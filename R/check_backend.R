#' @title Function \code{check_backend}
#' @description Checks that a given backend package.
#' If the package is not installed, the function quits in error.
#' @export
#' @param backend_package defaults to "fbseqCUDA". Other options include
#' "fbseqSerial", which does not use any parallel computing, and 
#' "fbseqOpenMP", which uses OpenMP for parallelization.
check_backend = function(backend_package = "fbseqCUDA"){
  if(!requireNamespace(backend_package, quietly = T))
    stop(paste0("function fbseq with specified backend requires package ", backend_package, "."))
}
