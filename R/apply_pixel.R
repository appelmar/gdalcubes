#' Apply a function over (multiband) pixels
#' 
#' This generic function applies a  function on pixels of a data cube, an R array, or other classes if implemented.
#' 
#' @param x object to be reduced 
#' @param ... additional arguments passed to method implementations
#' @return return value and type depend on the class of x
#' @seealso \code{\link{apply_pixel.cube}}
#' @seealso \code{\link{apply_pixel.array}} 
#' @export
apply_pixel <- function(x, ...) {
  UseMethod("apply_pixel")
}




#' Apply arithmetic expressions over all pixels of a data cube
#' 
#' Create a proxy data cube, which applies arithmetics expressions over all pixels of a data cube. Expressions may access band values by their name.
#'
#' @param x Source data cube
#' @param expr character vector with one or more arithmetic expressions (see Details)
#' @param names optional character vector with the same length as expr to specify band names for the output cube
#' @param ... not used
#' @return A proxy data cube object
#' @details gdalcubes uses the \href{https://github.com/ArashPartow/exprtk}{exprtk library} to evaluate expressions in C++, you can look at the library examples 
#' to see what kind of expressions you can execute.
#' @examples 
#'  L8_files <- list.files(system.file("L8NY18", package = "gdalcubes"), 
#'                         ".TIF", recursive = TRUE, full.names = TRUE)
#'  v = cube_view(l=388941.2, r=766552.4, b=4345299, t=4744931, 
#'          proj="EPSG:32618",
#'          nx = 497, ny=526, t0="2018-01", t1="2018-12", dt="P1M")
#'  L8.col = create_image_collection(L8_files, "L8_L1TP") 
#'  L8.cube = data_cube(L8.col, v) 
#'  L8.cube = select_bands(L8.cube, c("B04", "B05")) 
#'  L8.ndvi = apply_pixel(L8.cube, "(B05-B04)/(B05+B04)", "NDVI") 
#'  L8.ndvi
#'  L8.ndvi.median =  reduce_time(L8.ndvi, "median(NDVI)") 
#'  plot(L8.ndvi.median, key.pos=1, zlim=c(0,1))
#'  
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
apply_pixel.cube <- function(x, expr, names=NULL, ...) {
  stopifnot(is.cube(x))
  
  if (is.null(names)) {
    names <- paste("band", 1:length(expr), sep="")
  }
  
  x = libgdalcubes_create_apply_pixel_cube(x, expr, names)
  class(x) <- c("apply_pixel_cube", "cube", "xptr")
  return(x)
}



is.apply_pixel_cube  <- function(obj) {
  if(!("apply_pixel_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




