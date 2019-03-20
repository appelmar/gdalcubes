
#' Fill NA data cube pixels by time series interpolation
#' 
#' Create a proxy data cube, which fills NA pixels of a data cube by nearest neighbor or linear time series interpolation. 
#' 
#' @param cube source data cube
#' @param method interpolation method, can be "near" or "linear"
#' @return a proxy data cube object
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
fill_time <- function(cube, method="near") {
  stopifnot(is.cube(cube))
  
  x = libgdalcubes_create_fill_time_cube(cube, method)
  class(x) <- c("fill_time_cube", "cube", "xptr")
  return(x)
}



is.fill_time_cube  <- function(obj) {
  if(!("fill_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




