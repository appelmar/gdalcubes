#' Join bands of two identically shaped data cubes 
#' 
#' Create a proxy data cube, which joins the bands of two identically shaped data cubes. The resulting cube
#' will have bands from both input cubes.
#'
#' @param X first source data cube
#' @param Y second source data cube
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
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'                           srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v)
#' L8.cube.b04 = select_bands(raster_cube(L8.col, v), c("B04"))
#' L8.cube.b05 = select_bands(raster_cube(L8.col, v), c("B05"))
#' join_bands(L8.cube.b04,L8.cube.b05)
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details 
#' Names of bands will be taken from the input cubes. If both cubes, however, have bands with identical name, prefixes are added to all band names. Prefixes
#' are tried to derive from names of provided X and Y arguments (derived with \code{substitute}) or simply set to "X." and "Y.". 
#' @export 
join_bands <- function(X, Y) {
  stopifnot(is.cube(X))
  stopifnot(is.cube(Y))
  
  prefix_X = ""
  prefix_Y = ""
  if (anyDuplicated(c(names(X),names(Y))) > 0) {
    prefix_X = "X"
    prefix_Y = "Y"
    if (is.name(substitute(X))) {
      prefix_X = as.character(substitute(X))
    }
    if (is.name(substitute(Y))) {
      prefix_Y = as.character(substitute(Y))
    }
  }
  
  x = libgdalcubes_create_join_bands_cube(X, Y, prefix_X, prefix_Y)
  class(x) <- c("join_bands_cube", "cube", "xptr")
  return(x)
}



is.join_bands_cube  <- function(obj) {
  if(!("join_bands_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




