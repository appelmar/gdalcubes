#' Select bands of a data cube
#' 
#' Create a proxy data cube, which selects specific bands of a data cube. The resulting cube
#' will drop any other bands.
#'
#' @param cube source data cube
#' @param bands character vector with band names
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
#' \donttest{
#' plot(L8.rgb, rgb=3:1)
#' }
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @note For performance reasons, \code{select_bands} should always be called directly on a cube created with \code{\link{raster_cube}} and 
#' drop all unneded bands. This allows to reduce RasterIO and warp operations in GDAL.
#' @export
select_bands <- function(cube, bands) {
  stopifnot(is.cube(cube))
  
  x = libgdalcubes_create_select_bands_cube(cube, bands)
  class(x) <- c("select_bands_cube", "cube", "xptr")
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




