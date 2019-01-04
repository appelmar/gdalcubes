
#' Apply arithmetic expressions over all pixels of a data cube
#' 
#' Create a proxy data cube, which applies arithmetics expressions over all pixels of a data cube. Expressions may access band values by their name.
#'
#' @param cube Source data cube
#' @param expr character vector with one or more arithmetic expressions (see Details)
#' @param names optional character vector with the same length as expr to specify band names for the output cube
#' @return A proxy data cube object
#' @details gdalcubes uses the \href{https://github.com/ArashPartow/exprtk}{exprtk library} to evaluate expressions in C++, you can look at the library examples 
#' to see what kind of expressions you can execute.
#' @examples 
#' \dontrun{
#' gcbs_apply_pixel(XXXX, expr=c("(B08-B04)/(B08+B04)"))
#' }
#' @note This function returns a proxy object, i.e., it will not start any computations besides deriving the shape of the result.
#' @export
gcbs_apply_pixel <- function(cube, expr, names=NULL) {
  stopifnot(is.gcbs_cube(cube))
  
  if (is.null(names)) {
    names <- paste("band", 1:length(expr), sep="")
  }
  
  x = libgdalcubes_create_apply_pixel_cube(cube, expr, names)
  class(x) <- c("gcbs_apply_pixel_cube", "gcbs_cube", "xptr")
  return(x)
}



is.gcbs_apply_pixel_cube  <- function(obj) {
  if(!("gcbs_apply_pixel_cube" %in% class(obj))) {
    return(FALSE)
  }
  if (libgdalcubes_is_null(obj)) {
    warning("GDAL data cube proxy object is invalid")
    return(FALSE)
  }
  return(TRUE)
}




