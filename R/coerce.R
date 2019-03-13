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

  outnc = tempfile(fileext = ".nc")
  subdatasets = paste0("NETCDF:\"", outnc, "\":", names(from), sep="", collapse = NULL)
  
  write_ncdf(from, outnc)
  out = stars::read_stars(subdatasets)
  out = stars::st_set_dimensions(out, "x", point = FALSE)
  out = stars::st_set_dimensions(out, "y", point = FALSE)
  out = stars::st_set_dimensions(out, "time", point = FALSE, values=as.POSIXct(dimension_values(from, "S")$t, tz = "GMT"))
 
  return(out)
}