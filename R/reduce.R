


#' Reduce a data cube over the time dimension
#' 
#' Create a proxy data cube, which applies a reducer function over pixel time series of a data cube
#'
#' @param cube Source data cube
#' @param reducer Reducer function, currently "min", "max", "median", "mean", "count", "sd", "var", or "sum"
#' @return A proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
#' @examples 
#' \dontrun{
#' gcbs_reduce(XXXX, "min")}
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
gcbs_reduce <- function(cube, reducer=c("mean","median","min","max")) {
  stopifnot(is.gcbs_cube(cube))

  x = libgdalcubes_create_reduce_cube(cube, reducer)
  class(x) <- c("gcbs_reduce_cube", "gcbs_cube", "xptr")
  return(x)
}





is.gcbs_reduce_cube  <- function(obj) {
  if(!("gcbs_reduce_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




