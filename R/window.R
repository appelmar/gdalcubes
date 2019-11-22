#' Apply a moving window operation over time
#' 
#' This generic function applies a reducer function over a moving window over the time dimension of a data cube, an R array, or other classes if implemented.
#' @param x object to be reduced 
#' @param ... further arguments passed to specific implementations
#' @return value and type depend on the class of x
#' @seealso \code{\link{window_time.cube}} 
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
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-07"),
#'                           srs="EPSG:32618", nx = 400, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.nir = select_bands(L8.cube, c("B05"))
#' window_time(L8.nir, window = c(2,2), "min(B05)")  
#' window_time(L8.nir, kernel=c(-1,1), window=c(1,0))
#' 
#' \donttest{
#' plot(window_time(L8.nir, kernel=c(-1,1), window=c(1,0)), key.pos=1)
#' } 
#' @export
window_time <- function(x, ...) {
  UseMethod("window_time")
}


#' Apply a moving window function over the time dimension of a data cube
#' 
#' Create a proxy data cube, which applies one ore more moving window functions to selected bands over pixel time series of a data cube.
#' The fuction can either use a predefined agggregation function or apply a custom convolution kernel. 
#'
#' @param x source data cube
#' @param kernel numeric vector with elements of the kernel 
#' @param expr either a single string, or a vector of strings defining which reducers wlil be applied over which bands of the input cube
#' @param window integer vector with two elements defining the size of the window before and after a cell, the total size of the window is window[1] + 1 + window[2]
#' @param ... optional additional expressions (if expr is not a vector)
#' @return proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does).
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
#'                           bottom=4345299, top=4744931, t0="2018-01", t1="2018-07"),
#'                           srs="EPSG:32618", nx = 400, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.nir = select_bands(L8.cube, c("B05"))
#' L8.nir.min = window_time(L8.nir, window = c(2,2), "min(B05)")  
#' L8.nir.min
#' 
#' L8.nir.kernel = window_time(L8.nir, kernel=c(-1,1), window=c(1,0))  
#' L8.nir.kernel
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details 
#' The function either applies a kernel convolution (if the \code{kernel} argument is provided) or a general reducer function 
#' over moving temporal windows. In the former case, the kernel convolution will be applied over all bands of the input 
#' cube, i.e., the output cube will have the same number of bands as the input cubes. If a kernel is given and the \code{window} argument is missing, 
#' the window will be symmetric to the center pixel with the size of the provided kernel. For general reducer functions, the window argument must be provided and
#' several expressions can be used to create multiple bands in the output cube.
#' 
#' Notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median".
#'
#' 
#' @export
window_time.cube <- function(x, expr,  ..., kernel, window) {
  stopifnot(is.cube(x))
  
  
  if (!missing(kernel)) {
    if (!missing(expr)) {
      warning("argument expr will be ignored, applying kernel convolution")
    }
    if (length(list(...))> 0) {
      warning("additional arguments will be ignored, applying kernel convolution")
    }
    if (!missing(window)) {
      stopifnot(window[1] %% 1 == 0)
      stopifnot(window[2] %% 1 == 0)
      stopifnot(window[1] + 1 + window[2] == length(kernel))
    }
    else if (length(kernel) %% 2 == 0) {
      window =  c(length(kernel) / 2, length(kernel) / 2 - 1)
      warning(paste("length of kernel is even, using an asymmetric window (", window[1], ",", window[2], ")", sep=""))
    }
    else {
      # default
      window = rep((length(kernel) - 1) / 2 ,2)
    }
    x = libgdalcubes_create_window_time_cube_kernel(x, as.integer(window), as.double(kernel))
    class(x) <- c("window_time_cube", "cube", "xptr")
    return(x)
  }
  else {
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
    x = libgdalcubes_create_window_time_cube_reduce(x, as.integer(window), reducers, bands)
    class(x) <- c("window_time_cube", "cube", "xptr")
    return(x)
  }
  
 
  
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




