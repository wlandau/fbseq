#' @title Function \code{getDevice}
#' @description Get the index of the current CUDA-capable GPU.
#' @export
#' @return Integer index of the current CUDA-capable GPU, which is >= 0 and 
#' < number of devices.
getDevice = function(){
  check_fbseq_input("fbseqCUDA")
  fbseqCUDA::RgetDevice()
}

#' @title Function \code{getDeviceCount}
#' @description Get the number of CUDA-capable GPUs.
#' @export
#' @return Number of CUDA-capable GPUs.
getDeviceCount = function(){
  check_fbseq_input("fbseqCUDA")
  fbseqCUDA::RgetDeviceCount()
}

#' @title Function \code{setDevice}
#' @description Choose a CUDA-capable GPU to use for the MCMC.
#' @export
#' @return Integer index of the newly-set CUDA-capable GPU, which is >= 0 and 
#' < number of devices.
#' @param device Integer index of a CUDA-capable GPU. Must be >= 0 < number of GPUs.
setDevice = function(device){
  check_fbseq_input("fbseqCUDA")
  fbseqCUDA::RsetDevice(as.integer(device))
}
