#' Coerce gdalcubes object into a stars object
#' 
#' @param from object to coerce
#' @export
as_stars <- function(from) { 
  stopifnot(inherits(from, "cube"))
  if (!requireNamespace("stars", quietly = TRUE))
    stop("stars package not found, please install first") 

  fname = tempfile(fileext = ".nc")
  as_ncdf(from, fname)
  return(stars::read_stars(fname))
}