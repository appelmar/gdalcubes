#' gdalcubes: Earth Observation Data Cubes from Satellite Image Collections
#'
#' Processing collections of Earth observation images as on-demand multispectral, multitemporal raster data cubes. Users
#' define cubes by spatiotemporal extent, resolution, and spatial reference system and let 'gdalcubes' automatically apply cropping, reprojection, and 
#' resampling using the 'Geospatial Data Abstraction Library' ('GDAL'). Implemented functions on data cubes include reduction over space and time, 
#' applying arithmetic expressions on pixel band values, moving window aggregates over time, filtering by space, time, bands, and predicates on pixel values, 
#' exporting data cubes as 'netCDF' or 'GeoTIFF' files, and plotting.  The package implements lazy evaluation and 
#' multithreading. All computational parts are implemented in C++, linking to the 'GDAL', 'netCDF', 'CURL', and 'SQLite' libraries. 
#' See Appel and Pebesma (2019) <doi:10.3390/data4030092> for further details.
#'
#' @docType package
#' @name gdalcubes
#' @importFrom Rcpp sourceCpp
#' @useDynLib gdalcubes
#' 
#' @importFrom grDevices grey rainbow dev.off dev.size png col2rgb
#' @importFrom graphics axis box image.default layout lcm par plot rasterImage rect title legend lines
#' @importFrom stats quantile rnorm
#' @importFrom utils head download.file
#' @import RcppProgress jsonlite ncdf4 
#' 
#' 
#' 
NULL

