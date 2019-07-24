
#' Fill NA data cube pixels by simple time series interpolation
#' 
#' Create a proxy data cube, which fills NA pixels of a data cube by nearest neighbor or linear time series interpolation. 
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db")) 
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P3M", aggregation = "median")
#' L8.cube = raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.filled = fill_time(L8.rgb, "linear")
#' L8.filled
#' 
#' \donttest{
#' plot(L8.filled, rgb=3:1, zlim=c(5000,12000))
#' }
#' 
#' @details 
#' 
#' Please notice that completely empty (NA) time series will not be filled, i.e. the result cube might still contain NA values. 
#' 
#' @param cube source data cube
#' @param method interpolation method, can be "near" (nearest neighbor), "linear" (linear interpolation), "locf" (last observation carried forward), or "nocb" (next observation carried backward)
#' @return a proxy data cube object
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
fill_time <- function(cube, method="near") {
  stopifnot(is.cube(cube))
  
  x = libgdalcubes_create_fill_time_cube(cube, method)
  class(x) <- c("fill_time_cube", "cube", "xptr")
  return(x)
}



is.fill_time_cube  <- function(obj) {
  if(!("fill_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




