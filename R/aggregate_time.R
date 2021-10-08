
#' Aggregate data cube time series to lower temporal resolution
#' 
#' Create a proxy data cube, which applies an aggregation function over pixel time series to lower temporal resolution.
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
#' L8.two_monthly = aggregate_time(L8.rgb, "P6M", "min")
#' L8.two_monthly
#' 
#' \donttest{
#' plot(L8.two_monthly, rgb=3:1, zlim=c(5000,12000))
#' }
#' 
#' @details 
#' This function can be used to aggregate time series to lower resolution (similar to \link[raster]{aggregate}) or to regularize
#' a data cube with irregular (labeled) time axis. It is possible to change the unit of the temporal resolution (e.g. to create monthly composites from daily images).
#' The size of the cube may be expanded automatically if the original temporal extent is not divisible by the new temporal size of pixels.
#' 
#' 
#' @param cube source data cube
#' @param method aggregation method, one of "mean", "min", "max", "median", "count", "sum", "prod", "var", and "sd"
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @note This function is experimental and only supports method = "mean" as a proof of concept.
#' @export
aggregate_time <- function(cube, dt, method="mean") {
  stopifnot(is.cube(cube))
  
  x = libgdalcubes_create_aggregate_time_cube(cube, dt, method)
  class(x) <- c("aggregate_time_cube", "cube", "xptr")
  return(x)
}



is.aggregate_time_cube  <- function(obj) {
  if(!("aggregate_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




