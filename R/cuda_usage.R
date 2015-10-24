#' @title Function \code{check_heterosisCUDA}
#' @description Checks that package \code{heterosisCUDA} is installed. 
#' If the package is not installed, the function quits in error.
#' @export
check_heterosisCUDA = function(){
  if(!requireNamespace("heterosisCUDA", quietly = T))
    stop("this function requires a CUDA-capable NVIDIA GPU, a working installation of CUDA with the correct drivers, and the heterosisCUDA package.")
}

#' @title Function \code{getDevice}
#' @description Get the index of the current CUDA-capable GPU.
#' @export
#' @return Integer index of the current CUDA-capable GPU, which is >= 0 and 
#' < number of devices.
getDevice = function(){
  check_heterosisCUDA()
  heterosisCUDA::RgetDevice()
}

#' @title Function \code{getDeviceCount}
#' @description Get the number of CUDA-capable GPUs.
#' @export
#' @return Number of CUDA-capable GPUs.
getDeviceCount = function(){
  check_heterosisCUDA()
  heterosisCUDA::RgetDeviceCount()
}

#' @title Function \code{setDevice}
#' @description Choose a CUDA-capable GPU to use for the MCMC.
#' @export
#' @return Integer index of the newly-set CUDA-capable GPU, which is >= 0 and 
#' < number of devices.
#' @param device Integer index of a CUDA-capable GPU. Must be >= 0 < number of GPUs.
setDevice = function(device){
  check_heterosisCUDA()
  heterosisCUDA::RsetDevice(as.integer(device))
}
