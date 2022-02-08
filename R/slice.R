#' Extract a single time slice from a data cube
#' 
#' Create a proxy data cube, which extracts a time slice from a data cube defined by label (datetime string) or integer index.
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
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M", aggregation = "median")
#' L8.cube = raster_cube(L8.col, v, mask=image_mask("BQA", bits=4, values=16))
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.slice = slice_time(L8.rgb, "2018-03")
#' L8.slice
#' 
#' \donttest{
#' plot(L8.slice, rgb=3:1, zlim=c(5000,12000))
#' }
#' 
#' @details 
#' Either \code{datetime} or \code{it} must be non-NULL. If both arguments are provided, the integer index \code{it} is ignored.
#' 
#' @param cube source data cube
#' @param datetime character; datetime string of the time slice
#' @param it integer; index of the time slice (zero-based)
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @export
slice_time <- function(cube, datetime=NULL, it=NULL) {
  stopifnot(is.cube(cube))
  if (is.null(datetime) && is.null(it)) {
    stop("Missing required argument: either datetime or it must be provided")
  }
  if (!is.null(datetime) && !is.null(it)) {
    warning("Argument it will be ignored because datetime has been provided")
    it = NULL
  }
  
  if (!is.null(it)) {
    if (it %% 1 != 0) {
      stop("Invalid argument: it must be an integer number")
    }
    if (length(it) != 1) {
      stop("Invalid argument: length(it) is not equal to 1")
    }
    x = gc_create_slice_time_cube(cube, "", as.integer(it))
  }
  else {
    if (!is.character(datetime)) {
      stop("Invalid argument: datetime must be of type character")
    }
    if (length(datetime) != 1) {
      stop("Invalid argument: length(datetime) is not equal to 1")
    }
    x = gc_create_slice_time_cube(cube, datetime, 0)
  }
  class(x) <- c("slice_time_cube", "cube", "xptr")
  return(x)
}



is.slice_time_cube  <- function(obj) {
  if(!("slice_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}



#' Extract a single time series (spatial slice) from a data cube
#' 
#' Create a proxy data cube, which extracts a time series from a data cube defined by spatial coordinates or integer x and y indexes.
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
#' L8.ts = slice_space(L8.rgb, loc = c(5e05, 4400000))
#' L8.ts
#' 
#' \donttest{
#' plot(L8.ts, join.timeseries = TRUE)
#' }
#' 
#' @details 
#' Either \code{loc} or \code{i} must be non-NULL. If both arguments are provided, integer indexes \code{i} are ignored.
#' 
#' @param cube source data cube
#' @param loc numeric length-two vector; spatial coordinates (x, y) of the time series, expressed in the coordinate reference system of the source data cube
#' @param i integer length-2 vector; indexes (x,y) of the time slice (zero-based) 
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' 
#' @export
slice_space <- function(cube, loc=NULL, i=NULL) {
  stopifnot(is.cube(cube))
  if (is.null(loc) && is.null(i)) {
    stop("Missing required argument: either loc or i must be provided")
  }
  if (!is.null(loc) && !is.null(i)) {
    warning("Argument i will be ignored because loc has been provided")
    i = NULL
  }
  
  if (!is.null(i)) {
    if (i %% 1 != 0) {
      stop("Invalid argument: i must be exactly two integer number")
    }
    if (length(i) != 2) {
      stop("Invalid argument: length(i) is not equal to 2")
    }
    x = gc_create_slice_space_cube(cube, numeric(0), as.integer(i))
  }
  else {
    if (!is.numeric(loc)) {
      stop("Invalid argument: loc must be of type numeric")
    }
    if (length(loc) != 2) {
      stop("Invalid argument: length(loc) is not equal to 2")
    }
    x = gc_create_slice_space_cube(cube, loc, integer(0))
  }
  class(x) <- c("slice_space_cube", "cube", "xptr")
  return(x)
}



is.slice_space_cube  <- function(obj) {
  if(!("slice_space_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}






