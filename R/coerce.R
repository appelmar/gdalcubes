#' Coerce gdalcubes object into a stars object
#' 
#' The function materializes a data cube as a temporary netCDF file and loads the file 
#' with the stars package.
#' 
#' @param from data cube object to coerce
#' @return stars object
#' @examples 
#' \donttest{
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-04"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' as_stars(select_bands(raster_cube(L8.col, v), c("B04", "B05")))
#' }
#' @export
as_stars <- function(from) { 
  stopifnot(inherits(from, "cube"))
  if (!requireNamespace("stars", quietly = TRUE))
    stop("stars package not found, please install first") 

  outnc = tempfile(fileext = ".nc")
  #subdatasets = paste0("NETCDF:\"", outnc, "\":", names(from), sep="", collapse = NULL)
  
  write_ncdf(from, outnc)
  out = stars::read_ncdf(outnc)
  out = stars::st_set_dimensions(out, "x", point = FALSE)
  out = stars::st_set_dimensions(out, "y", point = FALSE)
  out = stars::st_set_dimensions(out, "time", point = FALSE, values=as.POSIXct(dimension_values(from, "S")$t, tz = "GMT"))
 
  return(out)
}