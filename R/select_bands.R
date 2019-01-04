#' Select bands of a data cube
#' 
#' Create a proxy data cube, which selects specific bands of a data cube. The resulting cube
#' will drop any other bands.
#'
#' @param cube Input data cube
#' @param bands character vector with band names
#' @return A proxy data cube object
#' @examples 
#' \dontrun{gcbs_select_bands()}
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @note For performance reasons, gcbs_select_bands should always be called directly on a image collection cube created with gcbs_cube and 
#' drop all unneded bands.This allows to reduce GDAL RasterIO and warp operations.
#' @export
gcbs_select_bands <- function(cube, bands) {
  stopifnot(is.gcbs_cube(cube))
  
  x = libgdalcubes_create_select_bands_cube(cube, bands)
  class(x) <- c("gcbs_select_bands_cube", "gcbs_cube", "xptr")
  return(x)
}




is.gcbs_select_bands_cube  <- function(obj) {
  if(!("gcbs_select_bands_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




