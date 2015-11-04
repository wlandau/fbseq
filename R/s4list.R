#' @title Function \code{s4list}
#' @description Turn an s4 object into a list
#' @export
#' @return a list of slots
#' @param obj an s4 object
s4list = function(obj){
  lst = lapply(slotNames(obj), function(x){slot(obj, x)})
  names(lst) = slotNames(obj)
  lst
}