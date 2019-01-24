


#' Reduce a data cube over the time dimension
#' 
#' Create a proxy data cube, which applies a reducer function over pixel time series of a data cube
#'
#' @param cube Source data cube
#' @param reducer Reducer function, currently "min", "max", "median", "mean", "count", "sd", "var", or "sum"
#' @return A proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = gcbs_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = gcbs_create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = gcbs_cube(L8.col, v) 
#'  L8.rgb = gcbs_select_bands(L8.cube, c("B02", "B03", "B04"))
#'  L8.rgb.median = gcbs_reduce(L8.rgb, "median")  
#'  L8.rgb.median
#'  plot(L8.rgb.median, rgb=3:1, zlim=c(4000,12000))
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @note This function is deprecated and will be replaced by the mor flexible gcbs_reduce_time.
#' @export
gcbs_reduce <- function(cube, reducer=c("mean","median","min","max")) {
  stopifnot(is.gcbs_cube(cube))

  x = libgdalcubes_create_reduce_cube(cube, reducer)
  class(x) <- c("gcbs_reduce_cube", "gcbs_cube", "xptr")
  return(x)
}





#' Reduce a data cube over the time dimension
#' 
#' Create a proxy data cube, which applies one ore more reducer functions over selected bands of pixel time series of a data cube
#'
#' @param cube Source data cube
#' @param expr Either a single string, or a vector of strings defining which reducers wlil be applied over which bands of the input cube
#' @param ... Optional additional expressions (if expr is not a vector)
#' @return A proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = gcbs_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = gcbs_create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = gcbs_cube(L8.col, v) 
#'  L8.rgb = gcbs_select_bands(L8.cube, c("B02", "B03", "B04"))
#'  L8.rgb.median = gcbs_reduce_time(L8.rgb, "median(B02)", "median(B03)", "median(B04)")  
#'  L8.rgb.median
#'  plot(L8.rgb.median, rgb=3:1, zlim=c(4000,12000))
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details Notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median", "var", "sd", "which_min", and "which_max".
#' @export
gcbs_reduce_time <- function(cube, expr, ...) {
  stopifnot(is.gcbs_cube(cube))
  stopifnot(is.character(expr))
  if (length(list(...))> 0) {
    stopifnot(all(sapply(list(...), is.character)))
    expr = c(expr, unlist(list(...)))
  }
  
  # parse expr to separate reducers and bands
  reducers = gsub("\\(.*\\)", "", expr)
  bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(reducers) == length(bands))
  x = libgdalcubes_create_reduce_time_cube(cube, reducers, bands)
  class(x) <- c("gcbs_reduce_time_cube", "gcbs_cube", "xptr")
  return(x)
}



#' Reduce a data cube over spatial (x,y or lat,lon) dimensions
#' 
#' Create a proxy data cube, which applies one ore more reducer functions over selected bands of spatial slices of a data cube
#'
#' @param cube Source data cube
#' @param expr Either a single string, or a vector of strings defining which reducers wlil be applied over which bands of the input cube
#' @param ... Optional additional expressions (if expr is not a vector)
#' @return A proxy data cube object
#' @note Implemented reducers will ignore any NAN values (as na.rm=TRUE does)
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"),
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = gcbs_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = gcbs_create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = gcbs_cube(L8.col, v) 
#'  L8.rgb = gcbs_select_bands(L8.cube, c("B02", "B03", "B04"))
#'  L8.rgb.median = gcbs_reduce_space(L8.rgb, "median(B02)", "median(B03)", "median(B04)")  
#'  L8.rgb.median
#'  plot(L8.rgb.median)
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @details Notice that expressions have a very simple format: the reducer is followed by the name of a band in parantheses. You cannot add
#' more complex functions or arguments.
#' 
#' Possible reducers currently are "min", "max", "sum", "prod", "count", "mean", "median", "var", "sd", "which_min", and "which_max".
#' @export
gcbs_reduce_space <- function(cube, expr, ...) {
  stopifnot(is.gcbs_cube(cube))
  stopifnot(is.character(expr))
  if (length(list(...))> 0) {
    stopifnot(all(sapply(list(...), is.character)))
    expr = c(expr, unlist(list(...)))
  }
  
  # parse expr to separate reducers and bands
  reducers = gsub("\\(.*\\)", "", expr)
  bands =  gsub("[\\(\\)]", "", regmatches(expr, gregexpr("\\(.*?\\)", expr)))
  stopifnot(length(reducers) == length(bands))
  x = libgdalcubes_create_reduce_space_cube(cube, reducers, bands)
  class(x) <- c("gcbs_reduce_space_cube", "gcbs_cube", "xptr")
  return(x)
}



is.gcbs_reduce_cube  <- function(obj) {
  if(!("gcbs_reduce_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.gcbs_reduce_time_cube  <- function(obj) {
  if(!("gcbs_reduce_time_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


is.gcbs_reduce_space_cube  <- function(obj) {
  if(!("gcbs_reduce_space_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}


