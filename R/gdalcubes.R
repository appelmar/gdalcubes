#' gdalcubes: Earth observation data cubes from GDAL image collections
#'
#' gdalcubes reads and processes collections of Earth Observation images as on-demand multispectral, multitemporal data cubes. Users
#' define cubes by spatiotemporal extent, resolution, and spatial reference system and let gdalcubes automatically apply cropping, reprojection, and 
#' resampling. Implemented functions on data cubes include reduction over space and time, applying arithmetic expressions on pixel band values,
#' moving window aggregates over time, filtering by space, time, bands, and predicates on pixel values, materializing data cubes as NetCDF files,
#' and plotting. Furthermore, user-defined R functions can be applied over chunks of data cubes. The package implements lazy evaluation and 
#' multithreading. All computational parts are implemented in C++, linking to the GDAL, NetCDF, CURL, and SQLite libraries.
#'
#' @docType package
#' @name gdalcubes
#' @importFrom Rcpp sourceCpp
#' @useDynLib gdalcubes
#' 
#' @importFrom grDevices grey
#' @importFrom graphics axis box image.default layout lcm par plot rasterImage rect title
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines
#' @importFrom stats quantile rnorm
#' @importFrom utils head 
#' @import RcppProgress jsonlite ncdf4 
#' 
#' 
#' 
NULL
