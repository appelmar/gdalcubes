#' Apply a moving window operation over time
#' 
#' This generic function applies a reducer function over a moving window over the time dimension of a data cube, an R array, or other classes if implemented.
#' @param x object to be reduced 
#' @param ... further arguments passed to specific implementations
#' @return Value and type depend on the class of x
#' @seealso \code{\link{window_time.cube}} 
#' @export
window_time <- function(x, ...) {
  UseMethod("window_time")
}


#' Apply a moving window function over the time dimension of a data cube
#' 
#' Create a proxy data cube, which applies one ore more moving window reducer functions over selected bands of pixel time series of a data cube
#'
#' @param x Source data cube
#' @param expr Either a single string, or a vector of strings defining which reducers wlil be applied over which bands of the input cube
#' @param window integer vector with two elements defining the size of the window before and after a cell, the total size ofof the window is window[1] + 1 + window[2]
#' @param ... Optional additional expressions (if expr is not a vector)
#' @return A proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
#' @examples 
#' L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                        ".TIF", recursive = TRUE, full.names = TRUE)
#' v = cube_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'               proj="EPSG:32618",
#'               nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#' L8.col = create_image_collection(L8_files, "L8_L1TP") 
#' L8.cube = data_cube(L8.col, v) 
#' L8.nir = select_bands(L8.cube, c("B08"))
#' L8.nir.min = window_time(L8.nir, c(2,2), "min(B02)")  
#' L8.nir.min
#' plot(L8.nir.min, zlim=c(4000,12000), key.pos=1, t=c(1,4,7))
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details Notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median".
#' @export
window_time.cube <- function(x, window, expr,  ...) {
  stopifnot(is.cube(x))
  stopifnot(is.character(expr))
  stopifnot(length(window) == 2)
  stopifnot(window[1] %% 1 == 0)
  stopifnot(window[2] %% 1 == 0)
  if (length(list(...))> 0) {
    stopifnot(all(sapply(list(...), is.character)))
    expr = c(expr, unlist(list(...)))
  }
  
  # parse expr to separate reducers and bands
  reducers = gsub("\\(.*\\)", "", expr)
  bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(reducers) == length(bands))
  x = libgdalcubes_create_window_time_cube(x, as.integer(window), reducers, bands)
  class(x) <- c("window_time_cube", "cube", "xptr")
  return(x)
}





is.window_time_cube  <- function(obj) {
  if(!("window_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




