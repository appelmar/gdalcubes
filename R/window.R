#' Apply a moving window operation over the time dimension of a data cube
#' 
#' Create a proxy data cube, which applies one ore more moving window functions to selected bands over pixel time series of a data cube.
#' The function can either apply a built-in aggregation function or apply a custom one-dimensional 
#' convolution kernel. 
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
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
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
#' Possible reducers include "min", "max", "sum", "prod", "count", "mean", and "median".
#'
#' 
#' @export
window_time <- function(x, expr,  ..., kernel, window) {
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
    x = gc_create_window_time_cube_kernel(x, as.integer(window), as.double(kernel))
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
    x = gc_create_window_time_cube_reduce(x, as.integer(window), reducers, bands)
    class(x) <- c("window_time_cube", "cube", "xptr")
    return(x)
  }
  
 
  
}





is.window_time_cube  <- function(obj) {
  if(!("window_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}



#' Apply a moving window (focal) operation or a convolution kernel over spatial dimensions of a data cube.
#' 
#' Create a proxy data cube, which applies a convolution kernel or aggregation functions over two-dimensional moving 
#' windows sliding over spatial slices of a data cube. The function can either execute one or more predefined aggregation functions or 
#' apply a custom convolution kernel. Among others, use cases include image processing (edge detection, noise reduction, etc.) and
#' enriching pixel values with local neighborhood properties (e.g. to use as predictor variables in ML models).
#'
#' @param x source data cube
#' @param kernel two dimensional kernel (matrix) applied as convolution (with odd number of rows and columns)
#' @param expr either a single string, or a vector of strings, defining which reducers will be applied over which bands of the input cube
#' @param window integer vector with two elements defining the size (number of pixels) of the window in y and x direction, the total size of the window is window[1] *  window[2]
#' @param keep_bands logical; if FALSE (the default), original data cube bands will be dropped. 
#' @param pad padding method applied to the borders; use NULL for no padding (NA), a numeric a fill value, or one of "REPLICATE", "REFLECT", "REFLECT_PIXEL"
#' @param ... optional additional expressions (if expr is not a vector)
#' @return proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as \code{na.rm = TRUE} does).
#' @note Calling this function consecutively many times may result in long computation times depending on chunk and window sizes due to the need to read adjacent data cube chunks.
#' @examples 
#' # create image collection from example Landsat data only 
#' # if not already done in other examples
#' if (!file.exists(file.path(tempdir(), "L8.db"))) {
#'   L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                          ".TIF", recursive = TRUE, full.names = TRUE)
#'   create_image_collection(L8_files, "L8_L1TP", file.path(tempdir(), "L8.db"), quiet = TRUE) 
#' }
#' 
#' L8.col = image_collection(file.path(tempdir(), "L8.db"))
#' v = cube_view(extent=list(left=388941.2, right=766552.4, 
#'                           bottom=4345299, top=4744931, t0="2018-04", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v, chunking = c(1,1000,1000)) 
#' L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#' L8.cube.mean5x5 = window_space(L8.cube, kernel = matrix(1/25, 5, 5))
#' L8.cube.mean5x5
#' 
#' \donttest{
#' plot(L8.cube.mean5x5, key.pos=1)
#' }
#' 
#' L8.cube.med_sd = window_space(L8.cube, "median(B04)" ,"sd(B04)", "median(B05)", "sd(B05)", 
#'                               window = c(5,5), keep_bands = TRUE)
#' L8.cube.med_sd
#' \donttest{
#' plot(L8.cube.med_sd, key.pos=1)
#' }
#' 
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details 
#' The function either applies a kernel convolution (if the \code{kernel} argument is provided) or one or more built-in reducer function 
#' over moving windows. 
#' 
#' In the former case, the kernel convolution will be applied over all bands of the input 
#' cube, i.e., the output cube will have the same number of bands as the input cubes.
#'  
#' To apply one or more aggregation functions over moving windows, the window argument must be provided as a vector with two integer sizes 
#' in the order y, x. Several string expressions can be provided to create multiple bands in the output cube. 
#' Notice that expressions have a very simple format: the reducer is followed by the name of a band in parentheses, e.g, "mean(band1)".
#' Possible reducers include "min", "max", "sum", "prod", "count", "mean", "median", "var", and "sd".
#' 
#' Padding methods "REPLICATE", "REFLECT", "REFLEX_PIXEL" are defined according to 
#' \url{https://openeo.org/documentation/1.0/processes.html#apply_kernel}.
#'
#' @export
window_space <- function(x, expr,  ..., kernel, window, keep_bands = FALSE, pad = NA) {
  stopifnot(is.cube(x))
  
  pad_fill = as.numeric(0)
  pad_mode = ""
  if (is.na(pad)) {
    pad_mode = ""
  }
  else if (is.numeric(pad)) {
    pad_fill = pad[1]
    pad_mode = "CONSTANT"
  }
  else if (is.character(pad)) {
    if (any(c("REPLICATE","replicate","edge") %in% pad)){
      pad_mode = "REPLICATE"
    }
    else if (any(c("REFLECT","reflect","symmetric") %in% pad)){
      pad_mode = "REFLECT"
    }
    else if (any(c("REFLECT_PIXEL","reflect_pixel") %in% pad)){
      pad_mode = "REFLECT_PIXEL"
    }
    else {
      warning("Unknown padding method (argument pad) provided: falling back to default method (no padding)")
      pad_mode = ""
    }
  }
  else {
    warning("Invalid padding method (argument pad) provided: falling back to default method (no padding)")
    pad_mode = ""
  }
  
  if (!missing(kernel)) {
    if (!missing(expr)) {
      warning("argument expr will be ignored, applying kernel convolution")
    }
    if (length(list(...))> 0) {
      warning("additional arguments will be ignored, applying kernel convolution")
    }
    if (!is.matrix(kernel)) {
      stop("Kernel must be provided as a matrix")
    }
    x = gc_create_window_space_cube_kernel(x, as.double(kernel), as.integer(nrow(kernel)), as.integer(ncol(kernel)), keep_bands, as.character(pad_mode), as.double(pad_fill))
  }
  else {
    stopifnot(is.character(expr))
    stopifnot(length(window) == 2)
    stopifnot(window[1] %% 1 == 0)
    stopifnot(window[2] %% 1 == 0)
    
    stopifnot(window[1] %% 2 == 1)
    stopifnot(window[2] %% 2 == 1)
    
    if (length(list(...))> 0) {
      stopifnot(all(sapply(list(...), is.character)))
      expr = c(expr, unlist(list(...)))
    }
    
    # parse expr to separate reducers and bands
    reducers = gsub("\\(.*\\)", "", expr)
    bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
    stopifnot(length(reducers) == length(bands))
    x = gc_create_window_space_cube_reduce(x, reducers, bands, as.integer(window[1]), as.integer(window[2]), keep_bands, as.character(pad_mode), as.double(pad_fill))
  }
  class(x) <- c("window_space_cube", "cube", "xptr")
  return(x)

}
  
  


is.window_space_cube  <- function(obj) {
  if(!("window_space_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (gc_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}





