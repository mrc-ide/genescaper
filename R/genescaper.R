#------------------------------------------------
#' @title insert title
#'
#' @description insert description
#'
#' @docType package
#' @name genescaper
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib genescaper, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("genescaper", libpath)
}
