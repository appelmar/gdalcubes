#' gdalcubes: Streaming GDAL Data Cubes to R
#'
#' The gdalcubes package constructs Earth Observation data cubes from many GDAL images / datasets and allows streaming data chunk-wise into R.
#'
#' @docType package
#' @name gdalcubes
#' @importFrom Rcpp sourceCpp
#' @useDynLib gdalcubes
#' 
#' @importFrom grDevices grey
#' @importFrom graphics axis box image.default layout lcm par plot rasterImage rect title
#' @importFrom stats quantile rnorm
#' @importFrom utils head 
#' @import RcppProgress jsonlite ncdf4 
#' 
#' 
NULL
