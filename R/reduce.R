#' Reduce multidimensional data over time
#' 
#' This generic function applies a reducer function over a data cube, an R array, or other classes if implemented.
#' @param x object to be reduced 
#' @param ... further arguments passed to specific implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{reduce_time.cube}} 
#' @seealso \code{\link{reduce_time.array}} 
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' reduce_time(raster_cube(L8.col, v) , "median(B02)", "median(B03)", "median(B04)")  
#'  
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' y <- reduce_time(x, function(v) {
#'   apply(v, 1, mean)
#' })
#'  
#' @export
reduce_time <- function(x, ...) {
  UseMethod("reduce_time")
}

#' Reduce multidimensional data over space
#' 
#' This generic function applies a reducer function over a data cube, an R array, or other classes if implemented.
#' @param x object to be reduced 
#' @param ... further arguments passed to specific implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{reduce_space.cube}} 
#' @seealso \code{\link{reduce_space.array}} 
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' reduce_space(raster_cube(L8.col, v) , "median(B02)")  
#' 
#' 
#' d <- c(4,16,32,32)
#' x <- array(rnorm(prod(d)), d)
#' y <- reduce_space(x, function(v) {
#'   apply(v, 1, mean)
#' })
#'  
#' @export
reduce_space <- function(x, ...) {
  UseMethod("reduce_space")
}

#' Reduce all bands of a data cube over the time dimension with a single reducer function
#' 
#' Create a proxy data cube, which applies a single reducer function over per-band pixel time series of a data cube
#'
#' @param cube source data cube
#' @param reducer reducer function, currently "min", "max", "median", "mean", "count", "sd", "var", or "sum"
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.rgb.median = reduce(L8.rgb, "median")  
#' L8.rgb.median
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @note This function is deprecated and will be replaced by the more flexible \code{\link{reduce_time}}.
#' @export
reduce <- function(cube, reducer=c("mean","median","min","max")) {
  stopifnot(is.cube(cube))

  x = libgdalcubes_create_reduce_cube(cube, reducer)
  class(x) <- c("reduce_cube", "cube", "xptr")
  return(x)
}


#' Reduce a data cube over the time dimension
#' 
#' Create a proxy data cube, which applies one or more reducer functions to selected bands over pixel time series of a data cube
#'
#' @param x source data cube
#' @param expr either a single string, or a vector of strings defining which reducers will be applied over which bands of the input cube
#' @param ... optional additional expressions (if \code{expr} is not a vector)
#' @return proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-06"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.rgb = select_bands(L8.cube, c("B02", "B03", "B04"))
#' L8.rgb.median = reduce_time(L8.rgb, "median(B02)", "median(B03)", "median(B04)")  
#' L8.rgb.median
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details Notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median", "var", "sd", "which_min", and "which_max".
#' @export
reduce_time.cube <- function(x, expr, ...) {
  stopifnot(is.cube(x))
  stopifnot(is.character(expr))
  if (length(list(...))> 0) {
    stopifnot(all(sapply(list(...), is.character)))
    expr = c(expr, unlist(list(...)))
  }
  
  # parse expr to separate reducers and bands
  reducers = gsub("\\(.*\\)", "", expr)
  bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(reducers) == length(bands))
  x = libgdalcubes_create_reduce_time_cube(x, reducers, bands)
  class(x) <- c("reduce_time_cube", "cube", "xptr")
  return(x)
}



#' Reduce a data cube over spatial (x,y or lat,lon) dimensions
#' 
#' Create a proxy data cube, which applies one or more reducer functions to selected bands over spatial slices of a data cube
#'
#' @param x source data cube
#' @param expr either a single string, or a vector of strings defining which reducers will be applied over which bands of the input cube
#' @param ... optional additional expressions (if \code{expr} is not a vector)
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
#'               bottom=4345299, top=4744931, t0="2018-01", t1="2018-12"),
#'               srs="EPSG:32618", nx = 497, ny=526, dt="P1M")
#' L8.cube = raster_cube(L8.col, v) 
#' L8.b02 = select_bands(L8.cube, c("B02"))
#' L8.b02.median = reduce_space(L8.b02, "median(B02)")  
#' L8.b02.median
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details Notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median", "var", "sd".
#' @export
reduce_space.cube <- function(x, expr, ...) {
  stopifnot(is.cube(x))
  stopifnot(is.character(expr))
  if (length(list(...))> 0) {
    stopifnot(all(sapply(list(...), is.character)))
    expr = c(expr, unlist(list(...)))
  }
  
  # parse expr to separate reducers and bands
  reducers = gsub("\\(.*\\)", "", expr)
  bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(reducers) == length(bands))
  x = libgdalcubes_create_reduce_space_cube(x, reducers, bands)
  class(x) <- c("reduce_space_cube", "cube", "xptr")
  return(x)
}



is.reduce_cube  <- function(obj) {
  if(!("reduce_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.reduce_time_cube  <- function(obj) {
  if(!("reduce_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.reduce_space_cube  <- function(obj) {
  if(!("reduce_space_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


