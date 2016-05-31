#' @title Function \code{backends}
#' @description Character vector of available MCMC backends for \code{fbseq}.
#' The elements are the R packages containing the backends, and the names
#' are the corresponding \code{backend} arguments to the \code{fbseq} function.
#' @export
#' @return Character vector of packages.
backends = function(){
  out = c("fbseqCUDA", "fbseqOpenMP")
  names(out) = gsub("fbseq", "", out)
  out
}
