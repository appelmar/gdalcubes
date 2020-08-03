#' Select time slices of a data cube
#' 
#' Create a proxy data cube, which selects specific time slices of a data cube. The time dimension of the resulting cube
#' will be irregular / labeled.
#'
#' @param cube source data cube
#' @param t character vector with date/time
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
#' L8.rgb = select_time(L8.rgb, c("2018-04", "2018-07"))
#' L8.rgb
#' \donttest{
#' plot(L8.rgb, rgb=3:1)
#' }
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
select_time <- function(cube, t) {
  stopifnot(is.cube(cube))
  
  x = libgdalcubes_create_select_time_cube(cube, t)
  class(x) <- c("select_time_cube", "cube", "xptr")
  return(x)
}




is.select_time_cube  <- function(obj) {
  if(!("select_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




