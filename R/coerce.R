#' Coerce gdalcubes object into a stars object
#' 
#' The function materializes a data cube as a temporary NetCDF file and loads the file 
#' with the stars package.
#' 
#' @param from data cube object to coerce
#' @return stars object
#' @export
as_stars <- function(from) { 
  stopifnot(inherits(from, "cube"))
  if (!requireNamespace("stars", quietly = TRUE))
    stop("stars package not found, please install first") 

  fname = tempfile(fileext = ".nc")
  write_ncdf(from, fname)
  return(stars::read_stars(fname))
}