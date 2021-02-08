#' Rename bands of a data cube
#' 
#' Create a proxy data cube, which renames specific bands of a data cube. 
#'
#' @param cube source data cube
#' @param ... named arguments with bands that will be renamed, see Details
#' @return proxy data cube object
#' @examples 
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
#'               bottom=4345299, top=4744931, t0="2018-04", t1="2018-07"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.rgb
#' L8.rgb = rename_bands(L8.cube, B02 = "blue", B03 = "green", B04 = "red")
#' L8.rgb
#' 
#' @details 
#' The result data cube always contains the same number of bands. No subsetting is 
#' done if only names for some of the bands are provided. In this case, only 
#' provided bands are renamed whereas other bands keep their original name. 
#' Variable arguments must be named by the old band name and the new names must be provided
#' as simple character values (see example). 
#' 
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
rename_bands <- function(cube, ...) {
  stopifnot(is.cube(cube))
  args = list(...)
  if (length(args) == 0) {
    return(cube)
  }
    
  arg_names = names(args)
  arg_values = unlist(args)
  
  stopifnot(is.character(arg_values))
  stopifnot(length(arg_names) == length(arg_values))
  
  x = libgdalcubes_create_rename_bands_cube(cube, arg_names, arg_values)
  class(x) <- c("rename_bands_cube", "cube", "xptr")
  return(x)
}




is.select_bands_cube  <- function(obj) {
  if(!("select_bands_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




