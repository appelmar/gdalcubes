
#' Spatial aggregation of data cubes
#' 
#' Create a proxy data cube, which applies an aggregation function to reduce the spatial resolution.
#' 
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", dx = 500, dy=500, dt="P3M", aggregation = "median")
#' L8.cube = raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.5km = aggregate_space(L8.rgb, 5000,5000, "mean")
#' L8.5km
#' 
#' \donttest{
#' plot(L8.5km, rgb=3:1, zlim=c(5000,12000))
#' }
#' 
#' @details 
#' This function reduces the spatial resolution of a data cube by applying an aggregation function to smaller blocks of pixels.   
#' 
#' The size of the cube may be expanded automatically in all directions if the original extent is not divisible by the new size of pixels.
#' 
#' Notice that if boundaries of the target cube do not align with the boundaries of the input cube (for example, if aggregating from 10m to 15m spatial resolution), pixels of
#' the input cube will contribute to the output pixel that contains its center coordinate. If the center coordinate is exactly on a boundary, the input pixel will contribute to 
#' the right / bottom pixel of the output cube. 
#' 
#' 
#' @param cube source data cube
#' @param dx numeric value; new spatial resolution in x direction
#' @param dy numeric value; new spatial resolution in y direction
#' @param method aggregation method, one of "mean", "min", "max", "median", "count", "sum", "prod", "var", and "sd"
#' @param fact simple integer factor defining how many cells (per axis) become aggregated to a single new cell, can be used instead of dx and dy
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @export
aggregate_space <- function(cube, dx, dy, method="mean", fact=NULL) {
  stopifnot(is.cube(cube))
  if (is.null(fact) && (missing(dx) || missing(dy))) {
    stop("Missing required argument: either dx and dy or fact must be provided")
  }
  if (!is.null(fact) && !missing(dx) && !missing(dy)) {
    warning("Argument fact will be ignored because dt has been provided")
    fact = NULL
  }
  
  if (!is.null(fact)) {
    if (fact %% 1 != 0) {
      stop("Invalid argument: fact must be an integer number > 1")
    }
    if (fact <= 1) {
      stop("Invalid argument: fact must be an integer number > 1")
    }
    x = gc_create_aggregate_space_cube(cube, as.double(0), as.double(0), method, as.integer(fact))
  }
  else {
    if (!is.numeric(dx)) {
      stop("Invalid argument: dx must be numeric")
    }
    if (!is.numeric(dy)) {
      stop("Invalid argument: dy must be numeric")
    }
    x = gc_create_aggregate_space_cube(cube, as.double(dx), as.double(dy), method, as.integer(0))
  }
  class(x) <- c("aggregate_space_cube", "cube", "xptr")
  return(x)
}



is.aggregate_space_cube  <- function(obj) {
  if(!("aggregate_space_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




