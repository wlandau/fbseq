#' @title Function \code{check_fbseq_input}
#' @description Checks input to the \code{fbseq} function.
#' @export
#' @return Corrected number of processes.
#' @param backend Backend for parallel computing. Options are
#' "CUDA" and "OpenMP".
#' @param processes number of CPU processes to fork for 
#' the additional chains. This argument is automatically reset to 1 for non-serial
#' backends because additional \code{parallel::mclapply} processes interfere
#' with other modes of parallelism.
#' For some other backends, chains will be distributed accross processes.
#' @param threads Number of threads for the OpenMP implementation.
check_fbseq_input = function(backend, processes = 1, threads = 1){
  backend_package = paste0("fbseq", backend)
  if(!(backend_package %in% backends())){
    stop("Illegal backend ", backend, ". Type backends() or ?backends for legal backends.")
  } else if(!requireNamespace(backend_package, quietly = T)){
    stop("function fbseq with ", backend, " backend requires package ", backend_package, ".")
  }

  if(threads > 1){
    if(backend == "CUDA"){
      warning("threads argument to the fbseq function is disregarded for the CUDA backend.")
    } else if (backend == "OpenMP" & !fbseqOpenMP::OpenMP_working()){
      warning("OpenMP is disabled, but more than 1 OpenMP thread requested.")
    }
  }

  if(processes > 1){
    if(backend == "CUDA"){
      processes = 1
      warning("For CUDA backend, processes must equal 1. Using 1 process.")
    } else if(threads > 1 & fbseqOpenMP::OpenMP_working()){
      processes = 1
      warning("For OpenMP backend and threads > 1, processes must equal 1. Using 1 process.")
    }
  }
  processes
}
