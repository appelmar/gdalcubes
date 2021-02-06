
#' Read a data cube from an existing netCDF file
#' 
#' Create a proxy data cube, which reads a data cube from a netCDF file that has been created using \code{\link{write_ncdf}}.
#' This function does not read cubes from arbitrary netCDF files and can be used e.g., to load intermediate results and/or 
#' plotting existing netCDF cubes on disk without doing the data cube creation from image collections.
#' 
#' @examples 
#' TODO
#'
#' 
#' @param path path to an existing netCDF file
#' @param chunking custom chunk sizes to read form the netCDF file; defaults to using chunk sizes from the netCDF file
#' @param auto_unpack logical; automatically apply offset and scale when reading data values
#' @return a proxy data cube object
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
ncdf_cube <- function(path, chunking = NULL, auto_unpack = TRUE) {
  stopifnot(file.exists(path))
  
  if (is.null(chunking)) {
    chunking = integer(0)
  }
  x = libgdalcubes_create_ncdf_cube(path, chunking, auto_unpack)
  class(x) <- c("ncdf_cube", "cube", "xptr")
  return(x)
}



is.ncdf_cube  <- function(obj) {
  if(!("ncdf_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




